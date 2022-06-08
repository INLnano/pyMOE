import matplotlib.pyplot as plt


from pyMOE import Aperture


def plot_aperture(aperture, scale=None, colorbar=True, only_plot=False, filename=None, **kwargs):
    """
    Plots the given aperture

    Args:
        aperture: Aperture object of the mask
        colorbar: True/False flag to plot colorbars
        only_plot: if True, only shows image without labels and axes
        filename: if provided, saves figure to filename
    """
    assert type(aperture) is Aperture, "aperture given is not an Aperture object"

    
    if scale is not None:
        scale_factor = scale
    else:
        scale_factor = 1 
        scale = ""

    fig = plt.figure()
    ax = plt.gca()
    
    plt.sca(ax)
    ax.set_aspect(1)
    pcm = plt.pcolormesh(aperture.x/scale_factor, aperture.y/scale_factor, aperture.aperture,)
    
    if not only_plot:
        plt.xlabel("x [%sm]"%scale)
        plt.ylabel("y [%sm]"%scale)
        if colorbar:
            fig.colorbar(pcm, ax=ax, label='z [a.u.]', shrink=0.6)

    else:
        plt.axis('off')
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    plt.subplots_adjust(wspace=0.3)
    if filename is not None:
        plt.savefig(filename)





def plot_field(aperture, which='both', scale=None, colorbar=True, only_plot=False, filename=None, **kwargs):
    """
    Plots the given aperture

    Args:
        aperture: Aperture object of the mask
        which: default "both", "amplitude" or "phase"
        colorbar: True/False flag to plot colorbars
        only_plot: if True, only shows image without labels and axes
        filename: if provided, saves figure to filename
    """
    assert type(aperture) is ApertureField, "aperture given is not an Aperture object"
    assert which in ["both", "amplitude", "phase"]
    
    if scale is not None:
        scale_factor = scale
    else:
        scale_factor = 1 
        scale = ""

    if which == "both":
        fig, axes = plt.subplots(1,2, sharey=True, sharex=True, )
        ax1 = axes[0]
        ax2 = axes[1]
    elif which == "amplitude":
        fig = plt.figure()
        ax1 = plt.gca()
    elif which == "phase":
        fig = plt.figure()
        ax2 = plt.gca()
    
    if which in ["both", "amplitude"]:
        plt.sca(ax1)
        ax1.set_aspect(1)
        pcm = plt.pcolormesh(aperture.x/scale_factor, aperture.y/scale_factor, aperture.amplitude,)
        
        if not only_plot:
            plt.xlabel("x [%sm]"%scale)
            plt.ylabel("y [%sm]"%scale)
            if colorbar:
                fig.colorbar(pcm, ax=ax1, label='Amplitude [a.u.]', shrink=0.6)
            else:
                plt.title("Amplitude")
        else:
            plt.axis('off')
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)

    if which in ["both", "phase"]:
        plt.sca(ax2)
        ax2.set_aspect(1)
        pcm = plt.pcolormesh(aperture.x/scale_factor, aperture.y/scale_factor, aperture.phase)
    
        if not only_plot:
            plt.xlabel("x [%sm]"%scale)
            plt.ylabel("y [%sm]"%scale)
            if colorbar:
                fig.colorbar(pcm, ax=ax2, label='Phase [rad]', shrink=0.6)
            else:
                plt.title("Phase")
        else:
            plt.axis('off')
            ax2.get_xaxis().set_visible(False)
            ax2.get_yaxis().set_visible(False) 
            
    plt.subplots_adjust(wspace=0.3)
    if filename is not None:
        plt.savefig(filename)
