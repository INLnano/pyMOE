"""
optimizer.py 
Module for optimization functionalities 
""" 

#using jaxlib==0.4.30

    
from scipy.optimize import leastsq, least_squares, differential_evolution, basinhopping, dual_annealing, minimize


from pyMOE.generate import create_empty_aperture, create_empty_aperture_from_aperture, circular_aperture
from pyMOE.field import Screen, create_empty_field_from_aperture, create_empty_field_from_field, generate_uniform_field, modulate_field
from pyMOE.plotting import plot_field 
from pyMOE.propagate import Bluestein

import numpy as np
import matplotlib.pyplot as plt 

#import jax.numpy.fft as sfft
#import jax 

from pyMOE.utils import progress_bar, Timer

    
def bluestein_czt_jax(x, f1, f2, fs, mout):
    """
    Bluestein from Hu et al. 2020 
    x =  field.field  * F  
    
    TO COMPLETE
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

    #apply calc to obtain the CZT 
    # Y-direction 
    fy1 = ystart + fsy/ 2
    fy2 = yend + fsy / 2
    gout = bluestein_czt_jax(gout, fy1, fy2, fsy, myout)

    # X-direction 
    fx1 = xstart + fsx / 2
    fx2 = xend + fsx / 2
    gout = bluestein_czt_jax(gout, fx1, fx2, fsx, mxout)
      
    #Eq 6, now gout is the result of the CZT
    gout = F0 * gout

    return gout
    
    

def Bluestein_jax (field, screen, wavelength, n=None):
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

    def update_screen_slice(screen, g1, z_i):
        # screen: (H, W, Z), g1: (H, W)
        H, W, Z = screen.shape
        # Make g1 shape (H, W, 1) to match the 3D screen
        g1_3d = jnp.expand_dims(g1, axis=-1)
        
        # Update screen at slice [:, :, z_i] using dynamic_update_slice
        return lax.dynamic_update_slice(screen, g1_3d, (0, 0, z_i))


    xlen,ylen,zlen = screen.XX.shape

    with Timer():
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
    import jax.numpy as np 

    z_field = 0 # the field source is assumed at z=0
    r = np.sqrt( (field.XX-x)**2 + (field.YY-y)**2 + (z_field-z)**2)
    r2 = r*r

    prop1 = np.exp(r*1.0j*k)/r2
    prop2 = z * k/(2*np.pi) *( 1/(r*k) - 1.0j)
    propE = field.field * prop1 * prop2

    # integrate over the input field and return field
    if simp2d==True or method=="simp2d": 
        Exyz = simpson2d(propE,field.x[0], field.x[-1], field.y[0], field.y[-1]) /(2*np.pi)
        #print("Simpson 2D method")
    elif method=="qmc":
        if sampler is None: 
            print("Please add a sample selection of points, e.g. Sobol.")
        Exyz = qmc_2d(propE,field.x[0], field.x[-1], field.y[0], field.y[-1], sampler) /(2*np.pi)
        #print("QMC 2D method")
    elif method=="trap":
        Exyz = np.trapezoid(np.trapezoid(propE, x=field.x, axis=0), x=field.y, axis=0)
    else: 
        Exyz = integrate.simpson(integrate.simpson(propE, field.x),field.y)/(2*np.pi) 

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
    
def scalable_angular_spectrum_method_jax(field, screen, z, wavelength, pad_factor, skip_final_phase=True, \
                                     crop=False):
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
    
    if z > z_limit:
        print("Zlimit is " +str(z_limit)+" but z is "+str(z) )
        

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
    else:
        psi_final = psi_p_final

    return psi_final


def SASM_jax(field, screen, wavelength, pad_factor = 2, crop = False):
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
    
    screen_field = jnp.zeros_like(screen.XX, dtype=jnp.complex64)
    
    def update_screen_slice(screen, g1, z_i):
        from jax import lax 
        # screen: (H, W, Z), g1: (H, W)
        H, W, Z = screen.shape
        g1 = jnp.transpose(g1)
        # Make g1 shape (H, W, 1) to match the 3D screen
        g1_3d = jnp.expand_dims(g1, axis=-1)
        
        # Update screen at slice [:, :, z_i] using dynamic_update_slice
        return lax.dynamic_update_slice(screen, g1_3d, (0, 0, z_i))
    
    xlen,ylen,zlen = screen.XX.shape
    
    g0  = scalable_angular_spectrum_method_jax(field, screen, screen.z[0], wavelength, pad_factor= pad_factor, crop=crop)
    
    print(np.shape(g0))
    
    x1 = jnp.linspace(jnp.min(screen.x), jnp.max(screen.x), jnp.shape(g0)[0])
    y1 = jnp.linspace(jnp.min(screen.y), jnp.max(screen.y), jnp.shape(g0)[1])
    #z = np.linspace(wavelength, 3500*micro, 100)
    z1 = screen.z
    
    #print(z1)

    newscreen = Screen(x1,y1,z1)

    with Timer():
        for z_i in range(zlen):
            z = screen.ZZ[:, :, z_i][-1][0]

            g1  = scalable_angular_spectrum_method_jax(field, screen, z, wavelength, pad_factor= pad_factor, crop = crop)
            
            #print(np.shape(g1))
            
            #newscreen.screen[:, :, z_i] = g1
            newscreen.screen = update_screen_slice(newscreen.screen, g1, z_i)

            
    return newscreen
    
    
def propagate(xar, aperture, screen, wavelength, mask_amp = None, circ_radius=None, input_field="uniform", E0=1, propagation_method="bluestein", pad_factor=2): 
    """
    This function is wrapper of propagate module functionalities, by default employing bluestein method
    to obtain the field at a screen. The field is propagated from aperture populated with xar phase values 
    These xar will change during optimization depending on a loss function being minimized.     
    
    TODO: Add ASM propagator 
    
    Args: 
        :xar:                2D real mesh grid of the phase values to populate the XY aperture phase values 
        :aperture:           Aperture object that will hold xar values, and modulate a field E 
        :screen:             Screen object with propagated field 
        :wavelength:         Wavelength of propagation 
        :circ_radius:        Radius of a circular aperture surrouding the xar values 
        :input_field:        Input field, by default "uniform field" 
        :E0:                 Input field amplitude, default = 1 
        :propagation_method: Propagation method from propagate module, default "bluestein" 
        '
    Returns: 
        :EXY:                Field XY at the screen at the last propagation z point 
    """
    
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
        
    if input_field=="uniform":
        # Generate a uniform field
        field = generate_uniform_field(field, E0=E0)
        
        
        #print(aperture1.aperture.shape, aperture3.aperture.shape)

        # Modulates the field 
        field = modulate_field(field, amplitude_mask=None, phase_mask=aperture3x, autograd=True )
        
    
        #moe.plotting.plot_field(field)
    
    elif input_field=="gaussian": 
        ##TO CORRECT 

        field = generate_uniform_field(field, E0=E0)
    
    #Propagate field to screen 
    if propagation_method=="nojax": 
        # Bluestein 
        EXYZ = Bluestein(field, screen, wavelength)
        
        #XY plane field -> In principle can be XYZ screen 
        EXY= (EXYZ.screen[:, :, -1])
        
    #Propagate field to screen 
    if propagation_method=="bluestein": 
        # Bluestein 
        EXYZ = Bluestein_jax(field, screen, wavelength)
        
        #XY plane field -> In principle can be XYZ screen 
        EXY= EXYZ[:, :, -1] 
        
    #Propagate field to screen 
    if propagation_method=="SASM": 
        # Bluestein 
        EXYZ = SASM_jax(field, screen, wavelength, pad_factor=pad_factor, crop =True)
        
        #XY plane field -> In principle can be XYZ screen 
        EXY= EXYZ.screen[:, :, -1] 
        
    if propagation_method=="RS": 
        # RS TO FIX 
        EXYZ = RS_integral_jax(field, screen, wavelength, parallel_computing = False, method="trap")
        EXY= EXYZ[:, :, -1] 
    
    ### TODO: OTHER PROPAGATION METHODS 

    # release memory 
    if aperture2x is not None: 
        del aperture2x 
    del aperture1x, aperture3x
    del field
    
    return EXY



def optimize(loss, x0, args1=None, optimizer_method="trf", ftol=1e-2, xtol=1e-8, gtol=1e-12, bounds =(-np.inf, np.inf), \
                niter =2, minimizer_kwargs=None, verbose=True, eps = None, learning_rate = 0.1, max_iters=100, max_nfev=1e6, \
                jax=True, *args, **kwargs): 
    """
    Optimize function using scipy optimizers, minizes the 'loss' function given as input 
    
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
    if jax==True: 
        import jax.numpy as np
        from jax import jacrev

        def cal_jac(x, *args): 
            return jacrev(lambda x: loss(x, *args))(x).ravel()
    else: 
        cal_jac = '3-point'
    
    if verbose==True: 
        v = 1
    elif verbose==False: 
        v = 0
    else: 
        v = 2
        
        
    if optimizer_method in ["adam", "rmsprop"]:
        import optax
    
        # Select optimizer
        optimizer = {"adam": optax.adam, "rmsprop": optax.rmsprop, "lf": optax.lbfgs}[optimizer_method](learning_rate)

        # Initialization
        x = np.array(x0)
        opt_state = optimizer.init(x)

        loss_fn = lambda x: loss(x, *args1)
        loss_and_grad_fn = (jax.value_and_grad(loss_fn))

        loss_history = []

        prev_loss = np.inf

        for i in range(max_iters):
            loss_val, grads = loss_and_grad_fn(x)

            # Convergence check
            if np.abs(prev_loss - loss_val) < ftol:
                if verbose:
                    print(f"Stopping at iteration {i}, loss change below tolerance ({ftol})")
                break

            updates, opt_state = optimizer.update(grads, opt_state, x)
            x = optax.apply_updates(x, updates)

            loss_history.append(loss_val)
            prev_loss = loss_val

        solution = x

    elif optimizer_method=='leastsq-lm': 
        solution = leastsq(loss, x0= x0,  args=args1 ,jac =cal_jac, ftol=ftol, xtol=xtol, gtol=gtol, verbose=v, max_nfev=max_nfev)
    elif (optimizer_method=='trf') or (optimizer_method=='dogbox'):  
        solution = least_squares(loss, x0=x0, args=args1 , jac =cal_jac,  ftol=ftol, xtol=xtol, gtol=gtol, method=optimizer_method, bounds=bounds, verbose=v, max_nfev=max_nfev, **kwargs)
    elif optimizer_method=='differential_evolution': 
        solution = differential_evolution(loss, x0=x0, jac =cal_jac, bounds = bounds, args =args1, verbose=v)
    elif optimizer_method=='dual_annealing': 
        solution = dual_annealing(loss, x0=x0, args =args1,jac =cal_jac,  bounds = bounds, verbose=v)
    elif optimizer_method=='basinhopping': 
        solution = basinhopping(loss, x0=x0,minimizer_kwargs=minimizer_kwargs, niter=max_iters)
    elif optimizer_method in ['Nelder-Mead', 'Powell','CG', 'BFGS', 'Newton-CG','L-BFGS-B', \
    'TNC', 'COBYLA', 'COBYQA', 'SLSQP', 'trust-constr', 'dogleg', 'trust-ncg', 'trust-exact','trust-krylov']:
        argas = args1
        solution = minimize(loss, x0, argas, jac =cal_jac, options={'gtol': ftol, 'disp': True, 'maxiter':50000},bounds=bounds)
    else: 
        print("Option optimizer_method= '"+str(optimizer_method)+"' not recognized")

    return solution


def generate_loss_function(metric): 
    """
    Generates a loss function for optimizer from a metric function (see metrics module), or addhoc 'metric' function 
    """
    return 