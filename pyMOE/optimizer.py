"""
optimizer.py 
Module for optimization functionalities 
""" 

#using jaxlib==0.4.30

    
from scipy.optimize import leastsq, least_squares, differential_evolution, basinhopping, dual_annealing, minimize

from pyMOE.generate import create_empty_aperture, create_empty_aperture_from_aperture, circular_aperture
from pyMOE.field import Field, Screen, create_empty_field_from_aperture, create_empty_field_from_field, generate_uniform_field, modulate_field
from pyMOE.plotting import plot_field 
from pyMOE.propagate import Bluestein
from pyMOE.ensemble import Ensemble 

import numpy as np
import matplotlib.pyplot as plt 

import jax.numpy.fft as sfft
import jax 

from pyMOE.utils import progress_bar, Timer

    
def bluestein_czt_jax(x, f1, f2, fs, mout):
    """
    Bluestein from Hu et al. 2020 
    """
    import jax.numpy.fft as sfft
    import jax.numpy as jnp
    
    m, n = x.shape #m = x, n = y 
    
    f11 = f1 + (mout * fs + f2 - f1) / (2 * mout)
    f22 = f2 + (mout * fs + f2 - f1) / (2 * mout)
    
    #Eq S15, complex starting point 
    a = jnp.exp(1j * 2 * jnp.pi * f11 / fs)
    #Eq S16, complex point spacing
    w = jnp.exp(-1j * 2 * jnp.pi / ((mout * fs) /((f22 - f11) )))
    
    
    #calculation of the size of the tile N
    mp = m + mout - 1
    N = jnp.power(2, jnp.ceil(jnp.log2(mp)).astype(jnp.int32))

    #defining the number of vals where to evaluate the czt, from -m +1 to m -1  (2m-1 size)
    h = jnp.arange(-m + 1, max(mout, m) )
    
    # W^(m^2 /2) 
    h_vals = w ** ((h ** 2) / 2)

    
    # W^(-m^2 /2) #all m, but cropped to N after 
    hin = 1 / h_vals[0:mp+1]
    
    ### fft ( W^(-m^2 /2) )
    ft = sfft.fft(hin , N) #if N>len(hin) zero pads the rest, increases the pad with mout 

    #b is A^(-m) * W^(m^2 /2) , positive m 
    b = a ** (-(jnp.arange(0,m))) * h_vals[m -1 : 2 * m-1]
    
    #apply the b in each tile * field 
    x_t = x * jnp.tile(b[:, jnp.newaxis], (1, n))
    
    ###fft ( field * A^(-m) * W^(m^2 /2)  )
    b_fft = sfft.fft(x_t, N, axis=0)

    ###convolution: ifft( fft[ field * A^(-m) * W^(m^2 /2)] * fft (W^(-m^2/2)) ) * W^(m^2/2)
    b_ifft = sfft.ifft(b_fft * jnp.tile(ft[:, jnp.newaxis], (1, n)), axis=0)
    bconv = b_ifft[m - 1:mp, :].T * jnp.tile(h_vals[m - 1:mp], (n, 1)) #positive m until mp, not 2m-1 
    
    #Eq S18, but in center of tiles 
    Mshift = (-m / 2) +  0.5 

    #Eq S19 , phase shift
    Pshift = jnp.exp(-1j * 2 * jnp.pi * (jnp.linspace(0, mout - 1, mout) * (f22 - f11) + mout *f11 ) * (Mshift) / (mout*fs))
    
    #tile over 
    Mshift_til = jnp.tile(Pshift, (n, 1))
    bout = bconv * Mshift_til 
    
    return bout
    
    
    
def scalar_bluestein(field, screen, d, wavelength):
    """
    Scalar diffraction computation using Bluestein's method.

    Args: 
            :field:         input Field
            :screen:        Observation Screen
            :d:             Observation distance z of the the Screen position 
            :wavelength:    Wavelength
            
    Returns:
            :screen:    Returns the screen populated with the result
    """
    k = 2 * np.pi / wavelength  

    xstart, xend = np.min(screen.x), np.max(screen.x)
    ystart, yend = np.min(screen.y), np.max(screen.y)
    
    x0, y0 = field.XX, field.YY
    x1, y1 = screen.XX[:,:,-1], screen.YY[:,:,-1]
    
    #Eqs 4 &5 
    F0 = np.exp(1j * k * d) / (1j * wavelength * d) * np.exp(1j * k / (2 * d) * (x1**2 + y1**2))
    F = np.exp(1j * k / (2 * d) * (x0**2 + y0**2))

    #arg of Eq 6 
    gout = field.field  * F
    
    #freqs for bluestein 
    fsx = wavelength * d / field.pixel_x
    fsy = wavelength * d / field.pixel_y
    
    mxout, myout = len(screen.x), len(screen.y)

    # X-direction 
    fx1 = xstart + fsx / 2
    fx2 = xend + fsx / 2
    gout = bluestein_czt_jax(gout, fx1, fx2, fsx, mxout)
      
    
    #apply calc to obtain the CZT 
    # Y-direction 
    fy1 = ystart + fsy/ 2
    fy2 = yend + fsy / 2
    gout = bluestein_czt_jax(gout, fy1, fy2, fsy, myout)

    #Eq 6, now gout is the result of the CZT
    gout = F0 * gout

    return gout
    
    

def Bluestein_jax (field, screen, wavelength, n=None, use_timer=True):
    """
    Implements the Bluestein propagation 
    
    Args: 
        
        :field:     input Field
        :screen:    Observation Screen
        :wavelength:    wavelength to consider
        :n:         refractive index of the propagation medium (default=1 for vacuum/air)
    Returns:
        :screen:    Returns the screen populated with the result
    """
    
    from jax import lax 
    import jax.numpy as jnp
    
    screen_field = jnp.zeros_like(screen.XX, dtype=jnp.complex64)  # or screen.screen if initialized

    xlen,ylen,zlen = screen.XX.shape

    with Timer(enabled=use_timer):
        for z_i in range(zlen):
            z = screen.ZZ[:, :, z_i][-1][0]

            g1  = scalar_bluestein(field, screen, z, wavelength)
            
            screen_field = update_screen_slice(screen_field, g1, z_i)

    return screen_field
    
    

def kernel_RS(field, k, x,y,z, simp2d=False, method=False, sampler=None):
    """
    Calculates the RS kernel integral from a field input aperture, assumed to be at z=0
    and returns the calculated E field
    
    Implements the Kernel in Mahajan 2011 part II eq 1-20 

    Args:
        :field:     input field
        :k:         Calculated wavenumber k=2pi/(wl*n)
        :x,y,z:     x, y, z coordinates of the screen point being evaluated
        :simp2d:    Defaults False, if True uses the simpson2d function
        :method     Choose between "simp2d", "qmc" or "trap", or defults to scipy simpson 
        :sampler:   Integer sampler for the QMC, e.g. sampler.integers(l_bounds=0, u_bounds=n, n=n) where n is side lenth 
    Returns:
        :E:         Calculated field
    """
    import jax.numpy as jnp 
    from scipy import integrate 

    z_field = 0 # the field source is assumed at z=0
    r = jnp.sqrt( (field.XX-x)**2 + (field.YY-y)**2 + (z_field-z)**2)
    r2 = r*r

    prop1 = jnp.exp(r*1.0j*k)/r2
    prop2 = z * k/(2*jnp.pi) *( 1/(r*k) - 1.0j)
    propE = field.field * prop1 * prop2

    # integrate over the input field and return field
    if simp2d==True or method=="simp2d": 
        Exyz = simpson2d(propE,field.x[0], field.x[-1], field.y[0], field.y[-1]) /(2*jnp.pi)
        #print("Simpson 2D method")
    elif method=="qmc":
        if sampler is None: 
            print("Please add a sample selection of points, e.g. Sobol.")
        Exyz = qmc_2d(propE,field.x[0], field.x[-1], field.y[0], field.y[-1], sampler) /(2*jnp.pi)
        #print("QMC 2D method")
    elif method=="trap":
        Exyz = jnp.trapezoid(jnp.trapezoid(propE, x=field.x, axis=0), x=field.y, axis=0)
    else: 
        Exyz = integrate.simpson(integrate.simpson(propE, field.y),field.x)/(2*jnp.pi) 

    return Exyz
    

def RS_integral_jax(field, screen, wavelength, n=None, simp2d=False, sampler=None, method=None, parallel_computing=False):
    """
    JAX-compatible Rayleigh-Sommerfeld integral propagation.
    
    Args:
        :field:         Input Field object
        :screen:        Observation Screen (assumed to hold .XX, .YY, .ZZ)
        :wavelength:    Wavelength
        :n:             Refractive index volume or scalar (optional)
    Returns:
        :screen_field:  Computed complex field at the screen
    """
    import jax.numpy as jnp
    from jax import vmap, jit
    from functools import partial

    k0 = 2 * jnp.pi / wavelength

    # Extract grid shapes
    xlen, ylen, zlen = screen.XX.shape

    # Initialize output screen
    screen_field = jnp.zeros((xlen, ylen, zlen), dtype=jnp.complex64)

    # Loop manually over indices (outer loop in Python OK; inner in JAX)
    for z_i in range(zlen):
        for y_i in range(ylen):
            for x_i in range(xlen):
                x = screen.XX[x_i, y_i, z_i]
                y = screen.YY[x_i, y_i, z_i]
                z = screen.ZZ[x_i, y_i, z_i]

                if n is not None:
                    k = 2 * jnp.pi * n[x_i, y_i, z_i] / wavelength
                else:
                    k = k0

                result = kernel_RS(field, k, x,y,z, simp2d, method, sampler, )

                screen_field = screen_field.at[x_i, y_i, z_i].set(result)

    return screen_field
    
    
def zero_pad_jax(arr, pad_factor=2):
    import jax.numpy as jnp 
    
    H, W = arr.shape
    H_new, W_new = int(pad_factor * H), int(pad_factor * W)
    out = jnp.zeros((H_new, W_new), dtype=arr.dtype)
    start_H = (H_new - H) // 2
    start_W = (W_new - W) // 2
    #out[start_H:start_H + H, start_W:start_W + W] = arr
    out = out.at[start_H:start_H + H, start_W:start_W + W].set(arr)
    return out

def crop_to_physical_size_jax(arr, dx_out, desired_size_metersx, dy_out, desired_size_metersy):
    """
    Center crop an array based on physical size and pixel spacing.
    """
    N_totalx, N_totaly = arr.shape
    crop_sizex, crop_sizey = int(desired_size_metersx / dx_out), int(desired_size_metersy / dy_out)
    startx = (N_totalx - crop_sizex) // 2
    starty = (N_totaly - crop_sizey) // 2
    #out = arr.at[start:start+crop_size, start:start+crop_size]
    return arr[startx:startx+crop_sizex, starty:starty+crop_sizey]
    
def resize_field_to_shape_jax(field, output_shape, method="linear"):
    """
    Resize a 2D complex array to the desired output shape by
    resizing real and imaginary parts separately.
    """
    import jax.numpy as jnp
    from jax import image as jimage
    
    amplitude = jnp.abs(field)
    phase = jnp.angle(field)

    amp_resized = jimage.resize(amplitude, output_shape, method=method)

    phase_resized = jimage.resize(phase, output_shape, method=method)

    return amp_resized * jnp.exp(1j * phase_resized)
    
    
    
def scalable_angular_spectrum_method_jax(field, screen, z, wavelength, pad_factor, skip_final_phase=True, \
                                     crop=True):
    """
    kernel based on Heintzmann et al. 2023  
    """
    import jax.numpy.fft as sfft
    import jax.numpy as jnp
    
    Nx, Ny = field.field.shape
    N = max([Nx,Ny]) 
    
    Lx = field.pixel_x * field.field.shape[0]  # Nx
    Ly = field.pixel_y * field.field.shape[1]  # Ny

    L = max([Lx, Ly]) 

    L_newx = pad_factor * Lx
    L_newy = pad_factor * Ly
    
    L_new = max([L_newx, L_newy])
    
    N_newx = int(pad_factor * Nx)
    N_newy = int(pad_factor * Ny)
    
    N_new = int(max([N_newx, N_newy]))
    
    psi_p = zero_pad_jax(field.field, pad_factor=pad_factor)
    
    
    z_limit = (- 4 * L * jnp.sqrt(8*L**2 / N**2 + wavelength**2) * jnp.sqrt(L**2 * 1 / (8 * L**2 + N**2 * wavelength**2))\
               / (wavelength * (-1+2 * jnp.sqrt(2) * jnp.sqrt(L**2 * 1 / (8 * L**2 + N**2 * wavelength**2)))))
    
    z_limit = (-4 * L_new * jnp.sqrt(8 * L_new**2 / N_new**2 + wavelength**2)
           * jnp.sqrt(L_new**2 / (8 * L_new**2 + N_new**2 * wavelength**2)) /
           (wavelength * (-1 + 2 * jnp.sqrt(2) * jnp.sqrt(L_new**2 / (8 * L_new**2 + N_new**2 * wavelength**2)))))
    
    #if z > z_limit:
    #    print("Zlimit is " +str(z_limit)+" but z is "+str(z) )
        

    k = 2 * jnp.pi / wavelength
    
    dfx = 1 / L_newx
    dfy = 1 / L_newy
    
    Lfx = N_newx * dfx
    Lfy = N_newy * dfy

    fx = sfft.fftfreq(N_newx, d=field.pixel_x)
    fy = sfft.fftfreq(N_newy, d=field.pixel_y)
    
    FX, FY = jnp.meshgrid(fx, fy, indexing='ij')

    x = sfft.ifftshift(np.linspace(-L_newx/2, L_newx/2, N_newx, endpoint=False))
    y = sfft.ifftshift(np.linspace(-L_newy/2, L_newy/2, N_newy, endpoint=False))
    X, Y = jnp.meshgrid(x, y, indexing='ij')

    cx = wavelength * FX
    cy = wavelength * FY
    tx = L_newx / (2 * z) + jnp.abs(cx)
    ty = L_newy / (2 * z) + jnp.abs(cy)

    W = ((cx**2 * (1 + tx**2) / tx**2 + cy**2 <= 1) &
         (cy**2 * (1 + ty**2) / ty**2 + cx**2 <= 1))

    H_AS = jnp.sqrt(jnp.clip(1 - (FX * wavelength)**2 - (FY * wavelength)**2, 0, None))
    H_Fr = 1 - (FX * wavelength)**2 / 2 - (FY * wavelength)**2 / 2
    delta_H = W * jnp.exp(1j * k * z * (H_AS - H_Fr))

    psi_fft = sfft.fft2(sfft.ifftshift(psi_p))
    psi_precomp = sfft.ifft2(psi_fft * delta_H)

    dqx = wavelength * z / L_newx
    dqy = wavelength * z / L_newy
    Qx = dqx * Nx * pad_factor
    Qy = dqy * Ny * pad_factor

    qx = jnp.fft.ifftshift(jnp.linspace(-Qx/2, Qx/2, N_newx, endpoint=False))
    qy = jnp.fft.ifftshift(jnp.linspace(-Qy/2, Qy/2, N_newy, endpoint=False))
    QX, QY = jnp.meshgrid(qx, qy, indexing='ij')

    H_1 = jnp.exp(1j * k / (2 * z) * (X**2 + Y**2))

    if skip_final_phase:
        psi_p_final = sfft.fftshift(sfft.fft2(H_1 * psi_precomp))
    else:
        H_2 = jnp.exp(1j * k * z) * jnp.exp(1j * k / (2 * z) * (QX**2 + QY**2))
        psi_p_final = sfft.fftshift(H_2 * sfft.fft2(H_1 * psi_precomp))

    #psi_final = zero_unpad(psi_p_final, field.field.shape, pad_factor=pad_factor)
    #psi_final = psi_p_final
    # After propagation:
    dx_out = wavelength * z / L_newx  # new pixel size
    dy_out = wavelength * z / L_newy  # new pixel size
    desired_output_sizex = jnp.max(screen.x) - jnp.min(screen.x) # e.g. 5 mm screen
    desired_output_sizey = jnp.max(screen.y) - jnp.min(screen.y) # e.g. 5 mm screen

    if crop ==True: 
        psi_final = crop_to_physical_size_jax(psi_p_final, dx_out, desired_output_sizex, dy_out, desired_output_sizey)
        psi_final = resize_field_to_shape_jax(psi_final, (N,N) )
    else:
        psi_final = psi_p_final

    return psi_final
    

def update_screen_slice(screen, g1, z_i):
    from jax import lax 
    import jax.numpy as jnp

    # screen: (H, W, Z), g1: (H, W)
    H, W, Z = screen.shape
    
    #g1 = jnp.transpose(g1)
    
    # Make g1 shape (H, W, 1) to match the 3D screen
    g1_3d = jnp.expand_dims(g1, axis=-1)
    
    # Update screen at slice [:, :, z_i] using dynamic_update_slice
    return lax.dynamic_update_slice(screen, g1_3d, (0, 0, z_i))
        


def SASM_jax(field, screen, wavelength, pad_factor = 2, crop =True, use_timer=True ):
    """
    Implements the Scalale Angular Spectrum Method propagation (Heintzmann et al. )
    
    Args: 
        :field:     input Field
        :screen:    Observation Screen
        :wavelength:    wavelength to consider
        :n:         refractive index of the propagation medium (default=1 for vacuum/air)
    Returns:
        :screen:    Returns the screen populated with the result
    """
    import jax.numpy as jnp
    from jax import lax
    
    #screen_field = jnp.zeros_like(screen.XX, dtype=jnp.complex64)
    
    xlen,ylen,zlen = screen.XX.shape
    
    g0  = scalable_angular_spectrum_method_jax(field, screen, screen.z[0], wavelength, pad_factor= pad_factor, crop=crop)
    
    #print(np.shape(g0))
    
    x1 = jnp.linspace(jnp.min(screen.x), jnp.max(screen.x), jnp.shape(g0)[0])
    y1 = jnp.linspace(jnp.min(screen.y), jnp.max(screen.y), jnp.shape(g0)[1])
    z1 = screen.z
    
    #print(z1)
    
    screen_field = jnp.zeros_like(g0, dtype=jnp.complex64)  # or screen.screen if initialized
    #newscreen = Screen(x1,y1,z1)

    with Timer(enabled=use_timer):
        for z_i in range(zlen):
            z = screen.ZZ[:, :, z_i][-1][0]

            g1  = scalable_angular_spectrum_method_jax(field, screen, z, wavelength, pad_factor= pad_factor, crop = crop)
            
            #print(np.shape(g1))
            
            #newscreen.screen[:, :, z_i] = g1
            #newscreen.screen = update_screen_slice(newscreen.screen, g1, z_i)
            screen_field = update_screen_slice(screen.screen, g1, z_i) 

            
    return screen_field
    
def zero_pad_jax(arr, pad_factor=2):
    import jax.numpy as jnp 
    
    H, W = arr.shape
    H_new, W_new = int(pad_factor * H), int(pad_factor * W)
    out = jnp.zeros((H_new, W_new), dtype=arr.dtype)
    start_H = (H_new - H) // 2
    start_W = (W_new - W) // 2
    #out[start_H:start_H + H, start_W:start_W + W] = arr
    out = out.at[start_H:start_H + H, start_W:start_W + W].set(arr)
    return out

def crop_to_physical_size_jax(arr, dx_out, desired_size_metersx, dy_out, desired_size_metersy):
    """
    Center crop an array based on physical size and pixel spacing.
    """
    N_totalx, N_totaly = arr.shape
    crop_sizex, crop_sizey = int(desired_size_metersx / dx_out), int(desired_size_metersy / dy_out)
    startx = (N_totalx - crop_sizex) // 2
    starty = (N_totaly - crop_sizey) // 2
    #out = arr.at[start:start+crop_size, start:start+crop_size]
    return arr[startx:startx+crop_sizex, starty:starty+crop_sizey]
    
def resize_field_to_shape_jax(field, output_shape, method="linear"):
    """
    Resize a 2D complex array to the desired output shape by
    resizing real and imaginary parts separately.
    """
    import jax.numpy as jnp
    from jax import image as jimage
    
    amplitude = jnp.abs(field)
    phase = jnp.angle(field)

    amp_resized = jimage.resize(amplitude, output_shape, method=method)

    phase_resized = jimage.resize(phase, output_shape, method=method)

    return amp_resized * jnp.exp(1j * phase_resized)
        

def ASM_kernel_jax(field, z, wavelength, input_extent, input_df, n = 1.0,  bandlimit = True, shift_yx = (0.0, 0.0),  kykx = (0.0, 0.0)):
    """
    Compute the ASM propagation kernel H(kx, ky)

    Args: 
        :field:        Input Field
        :z:            Distance to Screen plane 
        :wavelength:   Wavelength
        :input_extent: Spatial extent of the input field (might be already padded) 
        :n:            Refractive index of the propagation medium (default=1 for vacuum/air)
        :bandlimit:    Boolean, default False, if True enforces band limit filters akin Matushima & Shimobaba
        :shift_yx:     tuple (shift_y, shift_x) with shift at the output screen plane, default (0.,0.)  
        :kykx:         tuple (ky,kx) angular input direction, default = (0.,0.) 
    Returns:
        :H field kernel:        Returns the kernel 
        
    """
    import jax.numpy.fft as sfft
    import jax.numpy as jnp
    
    axes = (-2, -1) 
    
    Ny, Nx = field.Nx, field.Ny
    dy, dx = field.pixel_y, field.pixel_x
    fx, fy = field.fx, field.fy 
    
    fx_grid, fy_grid = field.FX, field.FY
    f_grid = jnp.stack((fx_grid, fy_grid), axis=-1)
    
    kykx_arr = jnp.asarray(kykx) / (2 * jnp.pi) 
    f_shifted = f_grid - kykx_arr
    f2_shifted = jnp.sum(f_shifted**2, axis=-1)
    
    phase_delay = jnp.sqrt(jnp.complex128(1 - (wavelength / n)**2 * f2_shifted))
        
    # shift in output plane
    shift_yx = jnp.asarray(shift_yx)
    out_shift = 2 * jnp.pi * jnp.sum(f_grid * shift_yx, axis=-1)
    
    phase = (2 * jnp.pi * (n / wavelength) * jnp.abs(z) * phase_delay) + out_shift
    
    propagator_field = jnp.where(z >= 0, jnp.exp(1j * phase), np.conj(jnp.exp(1j * phase)))
    
    if bandlimit:
        # Bandlimit (Matsushima & Shimobaba)
        z_arr = jnp.array(z)
        shift_yx_grid = shift_yx[jnp.newaxis, jnp.newaxis, :]
        input_extent_grid = input_extent[jnp.newaxis, jnp.newaxis, :]
        input_df_grid = input_df[jnp.newaxis, jnp.newaxis, :]

        k_limit_p = (
            ((shift_yx_grid + 1 / (2 * input_df_grid) ) ** -2 * z_arr**2 + 1) ** (-1 / 2)
        ) / wavelength * n
        
        k_limit_n = (
            ((shift_yx_grid - 1 / (2 * input_df_grid) ) ** -2 * z_arr**2 + 1) ** (-1 / 2)
        ) / wavelength * n

        # k0: Center of the bandlimit filter 
        k0 = (1 / 2) * (
            jnp.sign(shift_yx_grid + input_extent_grid) * k_limit_p
            + jnp.sign(shift_yx_grid - input_extent_grid) * k_limit_n
        )
        
        # k_width: Width of the bandlimit filter 
        k_width = (
            jnp.sign(shift_yx_grid + input_extent_grid) * k_limit_p
            - jnp.sign(shift_yx_grid - input_extent_grid) * k_limit_n
        )
        
        k_max = k_width / 2 # Half the width
        
        # H band limit filter 
        H_filter_yx = jnp.abs(f_grid - k0) <= k_max
        H_filter = H_filter_yx[..., 0] * H_filter_yx[..., 1]
        
        propagator_field = propagator_field * H_filter
        
    return sfft.ifftshift(propagator_field, axes=axes)


def czt_jax(x, m,  w, a, axis = -1):
    import jax.numpy.fft as sfft
    import jax.numpy as jnp
    
    n = x.shape[axis]
    n_czt = m + n - 1
    k = jnp.arange(n_czt)
    wk2 = w ** (k**2 / 2)
    Awk2 = a ** -k[:n] * wk2[:n]
    Fwk2 = sfft.fft(1 / jnp.hstack((wk2[n - 1 : 0 : -1], wk2[:m])), n_czt)
    wk2 = wk2[:m]

    x = jnp.moveaxis(x, axis, -1)
    y = sfft.ifft(sfft.fft(x * Awk2, n_czt, axis=-1) * Fwk2, axis=-1)
    y = y[..., n - 1 : n + m - 1] * wk2
    y = jnp.moveaxis(y, -1, axis)
    return y


def ASM_propagate_jax(field, screen, z, wavelength, pad_width, n = 1.0, mode = None, bl = True, shift = None, kykx = (0.0, 0.0)):
    """
    Angular Spectrum Method (ASM) propagation computation.

    Args: 
        :field:       Input Field
        :screen:      Observation Screen
        :z:           Distance to Screen plane 
        :wavelength:  Wavelength
        :pad_width:   Padding around area of interest (assumed constant all around), in nr. of pixels  
        :n:           refractive index of the propagation medium (default=1 for vacuum/air)
        :mode:        ASM mode, options: None (default) = convetional ASM; "czt" = Chirp Z-Transform; "BLAS" = Band-Limited ASM
        :bl:          Boolean, default False, if True enforces band limit filters akin Matushima & Shimobaba
        :shift:       tuple (shift_y, shift_x) with shift at the output screen plane, default None = calculates the shift from screen limits. If mode="czt" shift is not used. 
        :kykx:        tuple (ky,kx) angular input direction, default = (0.,0.) 
    Returns:
        :field:       Returns the calculated field
    """
    import jax.numpy.fft as sfft
    import jax.numpy as jnp
    
    # zero pad the field  --- asssumes symmetric padding!
    padding = ((pad_width, pad_width), (pad_width, pad_width))
    padded_vals = jnp.pad(field.field, padding, mode='constant', constant_values=0)

    #print(padded_vals.shape)

    paddedx = jnp.linspace(field.x.min() - pad_width*field.pixel_x, field.x.max() + pad_width*field.pixel_x, len(field.x) + 2*pad_width )
    paddedy = jnp.linspace(field.y.min() - pad_width*field.pixel_y, field.y.max() + pad_width*field.pixel_y, len(field.y) + 2*pad_width )

    #print(field.x, paddedx, field.x.shape, paddedx.shape)
    #print(field.y, paddedy, field.y.shape, paddedy.shape)
    
    padded_field = Field(paddedx, paddedy)
    
    padded_field.field = padded_vals
    
    axes = (-2, -1)  #define axes as last two dims for generalization, although with 2D does not make a difference
    spatial_shape = jnp.array(padded_field.shape)
    
    input_dx = np.array([field.pixel_x, field.pixel_y])
    mask_field_extent = input_dx * spatial_shape 
    input_df = 1.0 / mask_field_extent 
    

    # Handle shifts from screen coordinates or given by user 
    # NOTE, czt mode does not take shifts
    if (shift is None) & (mode != "czt"):
        ymin, ymax = jnp.min(screen.y), jnp.max(screen.y)
        xmin, xmax = jnp.min(screen.x), jnp.max(screen.x)

        shift_y = (ymax + ymin)/2
        shift_x = (xmax + xmin)/2
        
        shift_yx_for_kernel = (shift_x, shift_y) 
    
    elif (shift is None) & (mode == "czt"):
        shift_yx_for_kernel = (0.,0.) 
    elif (shift is not None):
        shift_yx_for_kernel = shift
    else: 
        shift_yx_for_kernel = (0.,0.)

    # Compute kernel H(kx, ky)
    kernel_H = ASM_kernel_jax(padded_field, z, wavelength, input_extent=mask_field_extent, \
                          input_df=input_df, shift_yx = shift_yx_for_kernel, bandlimit = bl, n = n, kykx=kykx )
    
    # APPLY KERNEL FORWARD TRANSFORM 
    field_transform = sfft.fftshift(sfft.fft2(sfft.ifftshift(padded_field.field, axes=axes), axes=axes) * kernel_H, axes=axes)

    # INVERSE TRANSFORM
    if mode=="czt":
        output_shape = jnp.array(screen.shape)
        ymin, ymax = screen.y.min(), screen.y.max()
        xmin, xmax = screen.x.min(), screen.x.max()
        
        output_dx = jnp.array([ (ymax - ymin) / (output_shape[0] - 1), (xmax - xmin) / (output_shape[1] - 1) ])
                
        # Scaling factor: alpha = output_dx / input_df
        alpha = output_dx / input_df 

        limits_min = [xmin, ymin]
        limits_max = [xmax, ymax]
        
        for d, axis in enumerate(axes):
            m = output_shape[d]
            
            # czt parameters, 
            # a: starting point on the circle (related to the min limit) w: angular step/ratio (related to the span/range)
            a_czt = jnp.exp(-1j * 2 * jnp.pi / mask_field_extent[d] * limits_min[d])
            w_czt = jnp.exp(1j * (2 * jnp.pi / mask_field_extent[d]) * (limits_max[d] - limits_min[d]) / (m - 1))

            # apply the czt
            field_czt = czt_jax(field_transform, m, w_czt, a_czt, axis=axis)
            
            # Ensure complex type just to be sure
            field_transform = field_czt.astype(field_transform.dtype)
            
            center = (m - 1) // 2 
            
            # phase compensation/modulation factor after the czt
            compensation = w_czt ** (-center * np.arange(m)) * (a_czt**center)
            
            field_transform = jnp.moveaxis(field_transform, axis, -1)
            field_transform = field_transform * compensation 
            field_transform = jnp.moveaxis(field_transform, -1, axis)
        
        # Apply the final scaling factor
        final_scaling = jnp.prod(1.0 / alpha)
        field_transform = field_transform * final_scaling
        
    elif (mode == "BLAS") & (bl==True):
        fx_grid, fy_grid = padded_field.FX, padded_field.FY
        f_grid = jnp.stack((fx_grid, fy_grid), axis=-1)
            
        output_shape = jnp.array(screen.shape)
        ymin, ymax = screen.y.min(), screen.y.max()
        xmin, xmax = screen.x.min(), screen.x.max()
        
        output_dx_y = (ymax - ymin) / (output_shape[0] - 1)
        output_dx_x = (xmax - xmin) / (output_shape[1] - 1)

        output_dx = jnp.array([output_dx_x, output_dx_y])
        
        # Scaling factor: alpha = output_dx / input_df (Eq 7)
        alpha = output_dx / input_df
        
        # Eq 9 of "Band-limited angular spectrum numerical propagation method
        # with selective scaling of observation window size and sample number"
        # (2012)
        wn = alpha * f_grid 
        
        # f = kernel for convolution, Eq 9 first term
        f = jnp.prod(jnp.exp(-1j * np.pi / alpha * wn**2), axis=-1)
        # B = modulated k-space field, Eq 9 second term
        B = field_transform * jnp.prod( (1 / alpha) * jnp.exp(1j * jnp.pi / alpha * wn**2), axis=-1)

        prefactor = jnp.prod( output_dx * jnp.exp(1j * jnp.pi / alpha * f_grid**2),axis=-1,)
        
        field_transform = prefactor * jax.scipy.signal.fftconvolve(B, f, mode="same", axes=axes)
    
        # crop to unpadded 
        y_slice = slice(int(padded_field.Ny/2) - int(screen.Ny/2),  int(padded_field.Ny/2) - int(screen.Ny/2) + screen.Ny)
        x_slice = slice(int(padded_field.Nx/2) - int(screen.Nx/2),  int(padded_field.Nx/2) - int(screen.Nx/2) + screen.Nx)
        field_transform = field_transform[x_slice, y_slice]

    else:   
        # just IFFT 
        propagated_field = sfft.fftshift(sfft.ifft2(sfft.ifftshift(field_transform, axes=axes), axes=axes), axes=axes)

        #fig = plt.figure() 
        #plt.imshow(np.abs(propagated_field))

        scaling_x = screen.pixel_x / padded_field.pixel_x
        scaling_y = screen.pixel_y / padded_field.pixel_y

        y_slice = slice(int(padded_field.Ny/2) - int(screen.Ny/2 *scaling_y),  \
                        int(padded_field.Ny/2) - int(screen.Ny/2 *scaling_y) + int(screen.Ny*scaling_y) )
        x_slice = slice(int(padded_field.Nx/2) - int(screen.Ny/2 *scaling_y), \
                        int(padded_field.Nx/2) - int(screen.Ny/2 *scaling_y) + int(screen.Nx*scaling_x) )
        field_transform = propagated_field[x_slice, y_slice]

        #fig = plt.figure() 
        #plt.imshow(np.abs(psi_final))
        
        if (screen.Nx != field.Nx) or (screen.Ny != field.Ny):
            field_transform = resize_field_to_shape_jax(field_transform, (screen.Nx,screen.Ny) )
        else: 
            field_transform = resize_field_to_shape_jax(field_transform, (field.Nx,field.Ny) )
        
    #print(field_transform.shape)
        
    return field_transform


def ASM_jax(field, screen, wavelength, pad, n = 1.0, mode = None, bl = True, shift = None, kykx = (0.,0.), use_timer=True ):
    """
    Implements the Angular Spectrum Method (ASM), with czt and bandlimit options.

    Args: 
        :field:       Input Field
        :screen:      Observation Screen
        :wavelength:  Wavelength
        :pad:         Padding around area of interest (assumed constant all around), in nr. of pixels  
        :n:           refractive index of the propagation medium (default=1 for vacuum/air)
        :mode:        ASM mode, options: None (default) = convetional ASM; "czt" = Chirp Z-Transform; "BLAS" = Band-Limited ASM
        :bl:          Boolean, default False, if True enforces band limit filters akin Matushima & Shimobaba
        :shift:       tuple (shift_y, shift_x) with shift at the output screen plane, default None = calculates the shift from screen limits. If mode="czt" shift is not used. 
        :kykx:        tuple (ky,kx) angular input direction, default = (0.,0.) 
    Returns:
        :screen:      Returns the screen populated with the result   
    """
    import jax.numpy as jnp
        
    xlen, ylen, zlen = screen.Nx, screen.Ny, screen.Nz

    screen_field = jnp.zeros_like(screen.XX, dtype=jnp.complex64)  # or screen.screen if initialized
    
    with Timer(enabled=use_timer):
        for z_i, z in enumerate(screen.z):
    
            g1 = ASM_propagate_jax(field, screen, z, wavelength, pad, n=n, mode=mode, bl=bl, shift=shift, kykx=kykx)
            
            if screen.Ny == 1:
                screen.screen = np.reshape(update_screen_slice(screen.screen, g1, z_i), (screen.Ny, screen.Nx))
            else:
                screen_field = update_screen_slice(screen.screen, g1, z_i)
            
            #progress_bar((z_i) / (zlen))
        #progress_bar(1)

    return screen_field    

    
     
   
def propagate_through_ensemble_jax(ensemble,  wavelength , corr=1, propagation_methods_array=None): 
    """
    Propagates though Ensemble object (various MOE surfaces)  
    
    Args:
        :ensemble:                   Ensemble object 
        :wavelength:                 Propagation wavelength
        :xar_plus_z:                 Array with optimization variables (and optionally z position as variable) 
        :propagation_methods_array:  Propagation method 
    Returns:
        :overall_field_at_screen:    Returns a screen concatenating both 
    """
    #initial aperture (the rest will be inside a loop)
    from jax import numpy as jnp 
    
    fieldi = ensemble.input_light_field 
    ensemble.field_array[0] = fieldi
    
    field0 = modulate_field(fieldi, amplitude_mask = ensemble.aperture_array_amp[0],\
                                      phase_mask=ensemble.aperture_array_phase[0], autograd=True)

    screen0 = ensemble.screen_array[0]
    
    if propagation_methods_array==None:
        propagation_methods_array= ["bluestein" for i in ensemble.screen_array]
    
    
    if propagation_methods_array[0]=="bluestein":
        EXYZ0 = Bluestein_jax(field0, screen0, wavelength)
    
        ensemble.field_array[1].field = EXYZ0[:,:, -1]/corr
        
        EXYZ = EXYZ0
        
        overall_field_at_screen = screen0
        overall_field_at_screen = EXYZ/jnp.max(EXYZ)

    elif propagation_methods_array[0]=="ASM": 
        EXYZ0 = ASM_jax(field0, screen0, wavelength, pad=int(len(field0.x) ))
        ensemble.field_array[1].field = EXYZ0[:,:, -1]/corr
        
        EXYZ = EXYZ0
        
        overall_field_at_screen = screen0
        overall_field_at_screen = EXYZ/jnp.max(EXYZ)

    elif propagation_methods_array[0]=="SASM": 
        EXYZ0 = SASM_jax(field0, screen0, wavelength, pad_factor=4, crop=True)
        ensemble.field_array[1].field = EXYZ0[:,:, -1]/jnp.max(EXYZ0[:,:, 0])
        
        EXYZ = EXYZ0
        
        overall_field_at_screen = screen0
        overall_field_at_screen = EXYZ/jnp.max(EXYZ)
        
    elif propagation_methods_array[0]=="rayleigh-sommerfeld":
        screen_XY = create_screen_XY(np.min(screen0.x), np.max(screen0.x), len(screen0.x), 
                                        np.min(screen0.y), np.max(screen0.y), len(screen0.y), 
                                        z=np.max(screen0.z) )

        EXYZ0 = RS_integral_jax(field0, screen_XY, wavelength, simp2d=True)
        ensemble.field_array[1].field = EXYZ0[:,:, -1]/jnp.max(EXYZ0[:,:, 0])
    
        overall_field_at_screen = screen0
        overall_field_at_screen.screen[:,:,-1] = EXYZ0.screen[:,:,0]/jnp.max(EXYZ0.screen[:,:, 0])
        
        #print("RS only calculated correctly the XY at the end Z, the other values of the screen object (e.g. YZ plane) have been repeated")
    
    EXYZ = EXYZ0
    overall_field_at_screen2 = EXYZ0
    
    #to the rest for all 
    for i in np.arange(1,len(ensemble.aperture_array_phase)): 
        #print(i)
        field = ensemble.field_array[i]
        
        aperture_amp = ensemble.aperture_array_amp[i]
        aperture_phase = ensemble.aperture_array_phase[i]
    
        field = modulate_field(field, amplitude_mask = aperture_amp, phase_mask=aperture_phase, autograd = True)
    
        ensemble.field_array[i].field = field
        
        screen = ensemble.screen_array[i]
        
        if propagation_methods_array[i]=="bluestein":
            EXYZ1 = Bluestein_jax(field, screen, wavelength)
            EXYZ = EXYZ1
            
            if i<len(ensemble.aperture_array_phase)-1 : 
                #ensemble.field_array[i+1].field = EXYZ0[:,:, -1]/jnp.max(EXYZ0[:,:, -1])
                ensemble.field_array[i+1].field = EXYZ[:, :, -1] / corr
            
            

            overall_field_at_screen = screen
            overall_field_at_screen = EXYZ/jnp.max(EXYZ)
            
        if propagation_methods_array[i]=="ASM":
            EXYZ1 = ASM_jax(field, screen, wavelength, pad=int(len(field0.x)) )
            EXYZ = EXYZ1
            
            if i<len(ensemble.aperture_array_phase)-1 : 
                #ensemble.field_array[i+1].field =EXYZ0[:,:, -1]/jnp.max(EXYZ0[:,:, -1])
                ensemble.field_array[i+1].field = EXYZ[:, :, -1] / corr
            
            

            overall_field_at_screen = screen
            overall_field_at_screen = EXYZ/jnp.max(EXYZ)

        if propagation_methods_array[i]=="SASM":
            EXYZ1 = SASM_jax(field, screen, wavelength, pad_factor=4, crop=True)
            EXYZ = EXYZ1
            
            if i<len(ensemble.aperture_array_phase)-1 : 
                #ensemble.field_array[i+1].field = EXYZ0[:,:, -1]/jnp.max(EXYZ0[:,:, -1])
                ensemble.field_array[i+1].field = EXYZ[:, :, -1] / corr
            
            EXYZ = EXYZ1

            overall_field_at_screen = screen
            overall_field_at_screen = EXYZ/jnp.max(EXYZ)
        
        elif propagation_methods_array[i]=="rayleigh-sommerfeld":
            if i==len(ensemble.aperture_array_phase)-1:
                screen_YZ = create_screen_YZ(np.min(screen.y), np.max(screen.y), len(screen.y), 
                                            np.min(screen.z), np.max(screen.z), len(screen.z), 
                                            #z=np.max(screen0.z) )
                                            x=0 )
                                            
                EXYZ0 = RS_integral_jax(field, screen_YZ, wavelength, simp2d=True)
                
                if i<len(ensemble.aperture_array_phase)-1 :
                    ensemble.field_array[i+1].field = EXYZ0.screen[:,:,-1]/jnp.max(EXYZ0.screen[:,:, -1])

                overall_field_at_screen = screen
                overall_field_at_screen.screen[int(len(screen0.x)/2),:,:] = EXYZ0.screen[:,0,:]
 
                EXYZ = screen

                
            else:  
                screen_XY = create_screen_XY(np.min(screen.x), np.max(screen0.x), len(screen0.x), 
                                                np.min(screen0.y), np.max(screen0.y), len(screen0.y), 
                                                z=np.max(screen0.z) )
                                                
                EXYZ0 = RS_integral_jax(field, screen_XY, wavelength, simp2d=True)

                if i<len(ensemble.aperture_array_phase)-1 :
                    ensemble.field_array[i+1].field = EXYZ0.screen[:,:,0]/jnp.max(EXYZ0.screen[:,:, -1])

                overall_field_at_screen = screen
                overall_field_at_screen.screen[:,:,-1] = EXYZ0.screen[:,:,0]

                #attributing this as screen to avoid dimension errors, however ONLY the XY is meaningful!
                #might be fixed 
                EXYZ = screen

        del field 

        #ensemble.screen_array[i] = EXYZ
        
        #print(EXYZ.screen.shape)        
        #overall_field_at_screen = screen0
        #overall_field_at_screen.screen = EXYZ.screen
            
        overall_field_at_screen2 = jnp.concatenate((overall_field_at_screen2, EXYZ), axis=2)
        
        
    return overall_field_at_screen2 

    

    
def propagate(xar, aperture, screen, wavelength, mask_amp=None, circ_radius=None, input_field="uniform", \
              input_field_amp=None, input_field_phase=None, E0=1, propagation_method="bluestein", \
              pad_factor=2, modedef = "czt", ensemble_mode=False, corr=1 , use_timer=True): 
    """
    This function is wrapper of propagate module functionalities, by default employing bluestein method
    to obtain the field at a 2D screen at the furthest z value (most common). 
    The field is propagated from the aperture populated with xar phase values.
    These xar will change during optimization depending on a loss function being minimized.     

    TODO: not restrict to the propagation to just the last 2D screen at z value, but use the 
          full calculated Screen(x,y,z) - although this could be also done at each loss function definition 
    
    Args: 
        :xar:                2D real mesh grid of the phase values to populate the XY aperture phase values (array of these for ensemble_mode)
        :aperture:           Aperture object that will hold xar values, and modulate a field E (array of these for ensemble_mode)
        :screen:             Screen object with propagated field (array of these for ensemble_mode)
        :wavelength:         Wavelength of propagation (array of these for ensemble_mode)
        :mask_amp:           2D real mesh grid of the phase values to populate the XY aperture amplitude values   
        :circ_radius:        Radius of a circular aperture surrouding the xar values 
        :input_field:        Input field, by default "uniform", if "custom" modulates field by 2D mesh grids amplitude & phase = input_field_amp & input_field_phase
        :input_field_amp:    2D real mesh grid of the amplitude values to modulate the field, default None 
        :input_field_phase:  2D real mesh grid of the phase values to modulate the field, default None  
        :E0:                 Input field amplitude, default = 1 
        :propagation_method: Propagation method from propagate module, default "bluestein" 
        :pad_factor:         Padding factor for SASM , default = 2 
        :modedef:            Mode for ASM = ["czt", "BLAS", None], default = "czt" 
        :ensemble_mode:      Ensemble mode is for ensemble propagation/optimization 
        :corr:               correction factor, default 1, if overflows are verified helps to play with it 
        
    Returns: 
        :EXY:                Field XY at the screen at the last propagation z point 
    """

    # To optimize the implementation -> it generates a new ensemble at each optimization step 
    if ensemble_mode==True:        
        aperture_array_amp = mask_amp #array of 2D matrices for the amp of each 
        aperture_array_phase = aperture #array of apertures, one for each surface in the path  
        screen_array = screen #array of screens, representing what comes after each surface 
        wavelength_array = wavelength #array, can  be single e.g. [lambda] 
        #xar is an array too, array of the 2D phases of the ensemble surfaces 
        
        prop_methods = propagation_method #array e-g. ["bluestein", "bluestein"]
        
        for ip, aperturex in enumerate(aperture_array_phase):
            # Create Apertur
            aperturein =  create_empty_aperture_from_aperture(aperturex)
            aperturein.aperture = xar[ip] 
            
            aperture_array_phase[ip] = aperturein  
            
        if input_field =="uniform": 
            field_input = create_empty_field_from_aperture(aperturex)
            field_input = generate_uniform_field(field_input, E0=1)
            input_light_field = field_input
            
        elif input_field=="custom": 
            # Modulates the field 
            field_input = create_empty_field_from_aperture(aperturex)
            field_input = generate_uniform_field(field_input, E0=1)

            if input_field_amp is not None: 
                aperture_amp = create_empty_aperture_from_aperture(aperturex)
                if input_field_amp.shape == aperture_amp.aperture.shape: 
                    aperture_amp.aperture = input_field_amp
                else: 
                    print("Shape of input field amplitude is NOT the same the aperture shape. ")
                
            if input_field_phase is not None:   
                aperture_phase = create_empty_aperture_from_aperture(aperturex)
                if input_field_phase.shape == aperture_phase.aperture.shape: 
                    aperture_phase.aperture = input_field_phase
                else: 
                    print("Shape of input field phase is NOT the same the aperture shape. ")

            #incoherent phase 
            elif input_field_phase=="random":
                aperture_phase = create_empty_aperture_from_aperture(aperturex)
                if input_field_phase.shape == aperture_phase.aperture.shape: 
                    aperture_phase.aperture = np.random.rand(np.flatten(aperture_phase.shape))*2*np.pi
                else: 
                    print("Shape of input field phase is NOT the same the aperture shape. ")
            

            input_light_field = modulate_field(field_input, amplitude_mask=aperture_amp, phase_mask=aperture_phase, autograd=True )
            
            del aperture_phase, aperture_amp
            #moe.plotting.plot_field(input_light_field)
        #print(screen_array)
    
        ensemble = Ensemble(aperture_array_amp, aperture_array_phase, screen_array,  wavelength_array, input_light_field)
        #print(ensemble.screen_array)
        
        EXY_array = [propagate_through_ensemble_jax(ensemble, wavelength, propagation_methods_array=prop_methods, corr=corr) for wavelength in wavelength_array]
        #print(ensemble.screen_array)
        
        EXY = EXY_array[-1][:, :, -1] 

    ###If ensemble mode is false, do simple 
    else:     
        # Create Aperture
        aperture1x = create_empty_aperture_from_aperture(aperture) 
        
        if circ_radius is not None:
            aperture2x = create_empty_aperture_from_aperture(aperture1x)
            mask_amp = circular_aperture(aperture2x, radius=circ_radius, center=(0,0))
        else: 
            aperture2x = None
        
        if mask_amp is not None: 
            mask_amp = mask_amp
            
        aperture3x =  create_empty_aperture_from_aperture(aperture1x)
        aperture3x.aperture = xar 
        
        ### Modulate field by the created aperture
        field = create_empty_field_from_aperture(aperture1x)
        
        del aperture1x

    if (input_field=="uniform") and (ensemble_mode==False):
        # Generate a uniform field
        field = generate_uniform_field(field, E0=E0)
        #print(aperture1.aperture.shape, aperture3.aperture.shape)

        # Modulates the field 
        field = modulate_field(field, amplitude_mask=None, phase_mask=aperture3x, autograd=True )
        #moe.plotting.plot_field(field)
        
        del aperture3x 
    
    elif input_field=="gaussian": 
        ##TO CORRECT 
        field = generate_uniform_field(field, E0=E0)
            
    else: 
        field = input_field
    
    #Propagate field to screen 
    if propagation_method=="nojax": 
        # Bluestein, no jax, just for debugging purposes 
        EXYZ = Bluestein(field, screen, wavelength, use_timer=use_timer)
        
        #XY plane field -> In principle can be XYZ screen 
        EXY= (EXYZ.screen[:, :, -1])
        
    #Propagate field to screen 
    if propagation_method=="bluestein": 
        # Bluestein jax
        EXYZ = Bluestein_jax(field, screen, wavelength, use_timer=use_timer)
        
        #XY plane field -> In principle can be XYZ screen 
        EXY= EXYZ[:, :, -1] 
        
    #Propagate field to screen 
    if propagation_method=="SASM": 
        # SASM
        EXYZ = SASM_jax(field, screen, wavelength, pad_factor=pad_factor, crop =True, use_timer=use_timer)
        
        #XY plane field -> In principle can be XYZ screen 
        EXY= EXYZ[:, :, -1] 
    
    
    #Propagate field to screen 
    if propagation_method=="ASM": 
        # ASM
        EXYZ = ASM_jax(field, screen, wavelength, pad=pad_factor, n = 1.0, mode = modedef, bl = True, shift = None, kykx = (0.,0.),  use_timer=use_timer)
        
        #XY plane field -> In principle can be XYZ screen 
        EXY= EXYZ[:, :, -1] 
        
        
    if propagation_method=="RS": 
        # RS TO FIX 
        EXYZ = RS_integral_jax(field, screen, wavelength,  use_timer=use_timer)
        EXY= EXYZ[:, :, -1] 
 
    #del aperture1x, aperture3x
    del field
    
    return EXY


def setup_optimizer_logger_batch(x_size, batch_size=1, log_dir="logs", name="optimize_log"):
    """
    Setup logger and memory structures for batch binary logging of optimizer variables
    
    Returns:
        :x_size:        size of x arrays (needed for reshape)
        :logger:        Python logger for text logs
        :logfile_txt:   path to text log file
        :batch_list:    list to accumulate batches
        :iter_counter:  dict for iteration counting 

    """

    import logging
    import os
    import numpy as np
    from datetime import datetime

    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    logfile_txt = os.path.join(log_dir, f"{name}_{timestamp}.txt")
    
    logfile_bin = os.path.join(log_dir, f"{name}_{timestamp}_x.npy")
    
    # Text logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(logfile_txt)
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    # Initialize batch list and counter
    batch_list = []
    iter_counter = {"k": 0}
    
    # Return info needed for optimizer
    return logger, logfile_txt, logfile_bin, batch_list, iter_counter, batch_size, fh


def optimize(loss, x0, args1=None, optimizer_method="trf", ftol=1e-2, xtol=1e-8, gtol=1e-12, bounds =(-np.inf, np.inf), \
                niter =2, minimizer_kwargs=None, verbose=True, eps = None, learning_rate = 0.1, max_iters=100, max_nfev=1e6, \
                jax=True, logger=None, logfile_bin=None, batch_list=None, batch_size=1, iter_counter=None, fh=None,\
                 *args, **kwargs):
    """
    Optimize function using scipy optimizers, minizes the 'loss' function given as input - to complete 
    
    Args: 
        :loss:               Loss function, to be minimized, having arguments 'x' (an array of optimization params) \
                              and possibly having arguments 'args1"  
        :x0:                 Initial array of values for 'x' that will be optimized  
        :args1:              Arguments of the loss function, given as list or array   
        :optimizer_method:   Optimizer method, local optimizers are "leastsq-lm", "trf", "dogbox",\
                              and global optimizers are "differential_evolution", "dual_annealing" or "basinhopping"  
        :tol:                Tolerance of the methods, defaults to 1e-2 
        :bounds:             Bounds for the optimizer methods, check if it is (-np.inf, np.inf) or an array for each param in the chosen optimizer_method
        
    Returns: 
        :solution:           The scipy solution object, to call the actual result do solution.x 
    """
    import numpy as onp
        
    if jax:
        import jax
        import jax.numpy as np
        from jax import jacrev

        def cal_jac(x, *args):
            return jacrev(lambda x: loss(x, *args))(x).ravel()
    else:
        import numpy as np
        cal_jac = '3-point'
    
    if verbose==True: 
        v = 1
    elif verbose==False: 
        v = 0
    else: 
        v = 2
        
    loss_history = []
    x_iter_history = []

    if iter_counter is None:
        iter_counter = {"k": 0}
    last_x = {"val": None}
    
    if batch_list is None:
        batch_list = []
    
    ###########################################33
    # Initiliazation of utils functions
    def loss_with_logging(x, *args):
        x_np = onp.asarray(x)
        k = iter_counter["k"]
        
        val = float(loss(x_np, *args, ))

        # Text log 
        x_str = onp.array2string(x_np, precision=6, suppress_small=True, threshold=1000)
        if logger is not None:
            logger.info("iter=%d | loss=%.6e | x=%s", k, val, x_str)

        # Append to batch
        batch_list.append(x_np.copy())

        # Write batch to .npy file if full
        if logfile_bin is not None and len(batch_list) >= batch_size:
            arr = onp.array(batch_list, dtype=onp.float64)
            with open(logfile_bin, "ab") as f:
                arr.tofile(f)
            batch_list.clear()

        # Save histories
        x_iter_history.append(x_np.copy())
        loss_history.append(val)

        iter_counter["k"] += 1
        last_x["val"] = x_np.copy()


        return loss(x, *args, )
        
    def scipy_like_result(x, fun, success, nfev):
        from scipy.optimize import OptimizeResult
        
        solution = OptimizeResult(x=x, fun=fun, success=success, nfev=nfev)
        
        return solution 
        
    ############################################## 
    # Choose the optimizer 
    
    # Stochastic optimizers 
    if optimizer_method == "cma-es":
        # METHOD CMA-ES https://pypi.org/project/cma/
        
        import cma
        
        opts = cma.CMAOptions() 
        opts['tolfun'] = ftol
        #opts['maxfevals'] = max_iters
        opts['maxiter'] = max_iters
        opts['bounds'] = [bounds[0], bounds[1]]

        es = cma.CMAEvolutionStrategy(x0, sigma0=kwargs.get("sigma0", 0.8), )

        while not es.stop():
            xs = es.ask()
            vals = [float(loss_with_logging(x, *args1)) for x in xs]
            es.tell(xs, vals)

            if verbose:
                es.disp()
                
            if es.result.iterations > max_iters or es.sigma < 1e-1:
                break

        res = es.result

        solution = scipy_like_result(x=res.xbest, fun=res.fbest, success=True, nfev=res.evaluations)
        
    elif optimizer_method == "spsa":
        # METHOD SPSA   https://www.jhuapl.edu/SPSA/ 
        # https://www.tutorialspoint.com/spsa-simultaneous-perturbation-stochastic-approximation-algorithm-using-python
        
        # SPSA parameters
        a = 0.7  ## these have been found empirically 
        c = 0.2

        def loss_spsa(x):
            return loss_with_logging(x, *args1)
            
        
        def SPSA(loss_function, theta, a, c, num_iterations, ftol=1e-6, xtol=None, window=100,):
            """
            A method for minimizing a loss function via simultaneous perturbation     stochastic approximation (SPSA).
            
            Parameters:
            loss_function (function): Loss function to be minimized.
            theta (numpy array): The initial values of the parameters.
            a (float): Size of the step when updating parameters.
            c (float): The noise parameter for randomly generating perturbations.
            num_iterations (int): Number of iterations in the algorithm.
            
            Returns:
            tuple: A tuple containing the best parameter values found and the 
            corresponding loss value.
            """

            import numpy as np

            theta = theta.copy()
            best_theta = theta.copy()
            best_loss = loss_function(theta)

            loss_history = [best_loss]
            stable_iters = 0

            for i in range(1, num_iterations + 1):

                delta = np.random.choice([-1, 1], size=theta.shape)

                loss_plus  = loss_function(theta + c * delta)
                loss_minus = loss_function(theta - c * delta)

                gradient = (loss_plus - loss_minus) / (2.0 * c * delta)

                theta_prev = theta.copy()
                theta = theta - a * gradient

                loss_curr = loss_function(theta)
                loss_history.append(loss_curr)

                # Track best
                if loss_curr < best_loss:
                    best_loss = loss_curr
                    best_theta = theta.copy()

                # ---- ftol: loss-based stopping (smoothed) ----
                if len(loss_history) >= window:
                    recent = np.mean(loss_history[-window:])
                    previous = np.mean(loss_history[-2*window:-window])

                    if abs(previous - recent) < ftol:
                        stable_iters += 1
                    else:
                        stable_iters = 0

                    if stable_iters >= window:
                        print(f"SPSA stopped by ftol at iteration {i}")
                        break

                # ---- xtol: parameter change ----
                if xtol is not None:
                    if np.linalg.norm(theta - theta_prev) < xtol:
                        print(f"SPSA stopped by xtol at iteration {i}")
                        break

            return best_theta, best_loss
            

        x_opt, lossopt = SPSA(loss_spsa, onp.asarray(x0, dtype=onp.float64), a, c, max_iters, ftol=ftol, xtol=xtol ) 
        
        solution = scipy_like_result(x=x_opt, fun=loss_spsa(x_opt), success=True, nfev=len(loss_history),)

    # JAX optimizers     
    elif optimizer_method in ["adam", "rmsprop", "lbfgsoptax", "adamw"]:
        def loss_jax(x, *args):
            return loss(x, *args,)  # must use jax.numpy internally
    
        # some fix for compatibility 
        import jax
        import jax.tree_util

        class _TreeShim:
            map = staticmethod(jax.tree_util.tree_map)

        jax.tree = _TreeShim
        
        import optax
    
        optimizer = {"adam": optax.adam, "rmsprop": optax.rmsprop, "lbfgsoptax": optax.lbfgs, "adamw": optax.adamw}[optimizer_method](learning_rate)
        loss_and_grad_fn = jax.value_and_grad(loss_jax)

        x = np.array(x0)
        opt_state = optimizer.init(x)

        prev_loss = onp.inf

        for i in range(max_iters):
            loss_val, grads = loss_and_grad_fn(x, *args1)

            x_np = onp.asarray(x)
            loss_f = float(loss_val)

            if logger is not None:
                logger.info("iter=%d | loss=%.6e | x=%s", i, loss_f, onp.array2string(x_np, precision=6))

            x_iter_history.append(x_np.copy())
            loss_history.append(loss_f)

            if onp.abs(prev_loss - loss_f) < ftol:
                if verbose:
                    print(f"Stopping at iteration {i}, loss change below tolerance")
                break

            updates, opt_state = optimizer.update(grads, opt_state, x)
            x = optax.apply_updates(x, updates)

            prev_loss = loss_f

        solution = scipy_like_result(x=onp.asarray(x), fun=loss_history[-1] if loss_history else None, success=True, nfev=len(loss_history),)
        
        solution.loss_history = loss_history
        solution.x_iter_history = onp.array(x_iter_history)
    
    # Scipy optimizers 
    #### Least squares optimizers 
    elif optimizer_method=='leastsq-lm': 
        solution = leastsq(loss_with_logging, x0= x0,  args=args1 ,jac =cal_jac, ftol=ftol, xtol=xtol, gtol=gtol, verbose=v, max_nfev=max_nfev)
    elif (optimizer_method=='trf') or (optimizer_method=='dogbox'):  
        solution = least_squares(loss_with_logging, x0=x0, args=args1 , jac =cal_jac,  ftol=ftol, xtol=xtol, gtol=gtol, method=optimizer_method, bounds=bounds, verbose=v, max_nfev=max_nfev, **kwargs)
    
    elif optimizer_method=='differential_evolution': 
        solution = differential_evolution(loss_with_logging, x0=x0, jac =cal_jac, bounds = bounds, args =args1, verbose=v)
    elif optimizer_method=='dual_annealing': 
        solution = dual_annealing(loss_with_logging, x0=x0, args =args1,jac =cal_jac,  bounds = bounds, verbose=v)
    elif optimizer_method=='basinhopping': 
        solution = basinhopping(loss_with_logging, x0=x0,minimizer_kwargs=minimizer_kwargs, niter=max_iters)

    elif optimizer_method in ['Nelder-Mead', 'Powell','CG', 'BFGS', 'Newton-CG','L-BFGS-B', \
    'TNC', 'COBYLA', 'COBYQA', 'SLSQP', 'trust-constr', 'dogleg', 'trust-ncg', 'trust-exact','trust-krylov']:
        def minimize_callback(xk):
            k = iter_counter["k"]
            x_np = onp.asarray(xk)
            val = float(loss(xk, *args1))

            if logger is not None:
                logger.info("iter=%d | loss=%.6e | x=%s", k, val, onp.array2string(x_np, precision=6, threshold=1000))

            x_iter_history.append(x_np.copy())
            loss_history.append(val)
            iter_counter["k"] += 1

            # Append to batch
            batch_list.append(x_np.copy())
            if logfile_bin is not None and len(batch_list) >= batch_size:
                arr = onp.array(batch_list, dtype=onp.float64)
                with open(logfile_bin, "ab") as f:
                    arr.tofile(f)
                batch_list.clear()
        
        argas = args1
        solution = minimize(loss_with_logging, x0, argas, jac =cal_jac, callback = minimize_callback, options={'ftol': ftol, 'disp': True, 'maxiter':50000},bounds=bounds)
      
    else: 
        print("Option optimizer_method= '"+str(optimizer_method)+"' not recognized")

    solution.loss_history = loss_history 
    solution.x_iter_history = onp.array(x_iter_history)
    
    ##TODO 
    # METHOD SPSA https://pypi.org/project/spsa/   https://www.jhuapl.edu/SPSA/ 
    
    ##TODO 
    # METHOD CMA-ES https://github.com/CMA-ES/pycma 
    
    if logger is not None: 
        logger.info("Optimization finished")
        fh.close()
        logger.removeHandler(fh)
    
    return solution

def generate_loss_function(metric): 
    """
    Generates a loss function for optimizer from a metric function (see metrics module), or addhoc 'metric' function 
    """
    return 