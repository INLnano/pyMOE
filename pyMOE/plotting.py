"""
plotting.py 
Module for plotting apertures ans save them as image files 


"""

import matplotlib.pyplot as plt
import cv2 
import numpy as np

from pyMOE import Aperture
from pyMOE import Field
from pyMOE import Screen
from pyMOE import ApertureField


def save_mask_plot(maskcir, xsize, ysize, filename):
    plt.ioff()
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    fig.set_dpi(1)

    plt.axis('equal')
    plt.axis('off')
    
    maskcir = np.flip(maskcir,0)
    plt.imshow(maskcir, vmin=0, vmax=1,extent =[0,xsize,0,ysize], cmap=plt.get_cmap("Greys"))
    
    plt.tight_layout(pad=0, w_pad=0, h_pad=0) 
    fig.tight_layout(w_pad=0, h_pad=0, pad =0)
    
    plt.gca().set_axis_off()
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0) 
    plt.margins(0,0)
    
    fig.set_size_inches((xsize), (ysize))
    fig.savefig("temp.png", bbox_inches=0,pad_inches = 0, dpi=1)
    plt.close()
    
    #remove any white padding
    im = cv2.imread("temp.png")
    cv2.imwrite(filename, im)
    
    plt.ion()
    

def plot_aperture(aperture, scale=None, colorbar=True, only_plot=False, filename=None, **kwargs):
    """
    Plots the given aperture

    Args:
        :aperture:  Aperture object of the mask
        :colorbar:  True/False flag to plot colorbars
        :only_plot: if True, only shows image without labels and axes
        :filename:  if provided, saves figure to filename 
        ###available output file extensions are same as opencv https://docs.opencv.org/3.4/d4/da8/group__imgcodecs.html 
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
            fig.colorbar(pcm, ax=ax, shrink=0.6)

    else:
        plt.axis('off')
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    plt.subplots_adjust(wspace=0.3)
    if filename is not None:
        plt.savefig("temp.png", bbox_inches='tight', pad_inches = 0)
        plt.close(fig1)
        img = cv2.imread("temp.png")
        cv2.imwrite(filename, img)
    plt.show()





def plot_field(field, which='both', scale=None, colorbar=True, only_plot=False, filename=None, **kwargs):
    """
    Plots the given field

    Args:
        :field:  Field object of the mask
        :which:     default "both", "amplitude" or "phase"
        :colorbar:  True/False flag to plot colorbars
        :only_plot: if True, only shows image without labels and axes
        :filename:  if provided, saves figure to filename
        ###available output file extensions are same as opencv https://docs.opencv.org/3.4/d4/da8/group__imgcodecs.html
    """
    assert type(field) is Field, "field given is not an Field object"
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
        pcm = plt.pcolormesh(field.x/scale_factor, field.y/scale_factor, field.amplitude,)
        
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
        pcm = plt.pcolormesh(field.x/scale_factor, field.y/scale_factor, field.phase)
    
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
        plt.savefig("temp.png", bbox_inches='tight', pad_inches = 0)
        plt.close(fig1)
        # img = cv2.imread("temp.png")
        # cv2.imwrite(filename, img)





def plot_screen_XY(screen, which='both', scale=None, colorbar=True, only_plot=False, filename=None, **kwargs):
    """
    Plots the given screen

    Args:
        :screen:  Screen object of the mask
        :which:     default "both", "amplitude" or "phase"
        :colorbar:  True/False flag to plot colorbars
        :only_plot: if True, only shows image without labels and axes
        :filename:  if provided, saves figure to filename
        ###available output file extensions are same as opencv https://docs.opencv.org/3.4/d4/da8/group__imgcodecs.html
    """
    assert type(screen) is Screen, "screen given is not an Screen object"
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
        pcm = plt.pcolormesh(screen.x/scale_factor, screen.y/scale_factor, screen.amplitude.reshape((len(screen.x),len(screen.y))))
        
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
        pcm = plt.pcolormesh(screen.x/scale_factor, screen.y/scale_factor, screen.phase.reshape((len(screen.x),len(screen.y))))
    
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
        plt.savefig("temp.png", bbox_inches='tight', pad_inches = 0)
        plt.close(fig1)
        # img = cv2.imread("temp.png")
        # cv2.imwrite(filename, img)


def plot_screen_YZ(screen, which='both', scale=None, colorbar=True, only_plot=False, filename=None, **kwargs):
    """
    Plots the given screen

    Args:
        :screen:  Screen object of the mask
        :which:     default "both", "amplitude" or "phase"
        :colorbar:  True/False flag to plot colorbars
        :only_plot: if True, only shows image without labels and axes
        :filename:  if provided, saves figure to filename
        ###available output file extensions are same as opencv https://docs.opencv.org/3.4/d4/da8/group__imgcodecs.html
    """
    assert type(screen) is Screen, "screen given is not an Screen object"
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
        # ax1.set_aspect(1)
        pcm = plt.pcolormesh(screen.z/scale_factor, screen.y/scale_factor, screen.amplitude[:,0])
        
        if not only_plot:
            plt.xlabel("z [%sm]"%scale)
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
        # ax2.set_aspect(1)
        pcm = plt.pcolormesh(screen.z/scale_factor, screen.y/scale_factor, screen.phase[:,0])
    
        if not only_plot:
            plt.xlabel("z [%sm]"%scale)
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
        plt.savefig("temp.png", bbox_inches='tight', pad_inches = 0)
        plt.close(fig1)
        # img = cv2.imread("temp.png")
        # cv2.imwrite(filename, img)



def plot_screen_ZZ(screen, which='both', scale=None, only_plot=False, filename=None, **kwargs):
    """
    Plots the given screen

    Args:
        :screen:  Screen object of the mask
        :which:     default "both", "amplitude" or "phase"
        :colorbar:  True/False flag to plot colorbars
        :only_plot: if True, only shows image without labels and axes
        :filename:  if provided, saves figure to filename
        ###available output file extensions are same as opencv https://docs.opencv.org/3.4/d4/da8/group__imgcodecs.html
    """
    assert type(screen) is Screen, "screen given is not an Screen object"
    assert which in ["both", "amplitude", "phase"]
    
    if scale is not None:
        scale_factor = scale
    else:
        scale_factor = 1 
        scale = ""

    if which == "both":
        fig, axes = plt.subplots(1,2, sharey=False, sharex=True, )
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
        # ax1.set_aspect(1)
        pcm = plt.plot(screen.z/scale_factor, screen.amplitude[0,0,:],)
        
        if not only_plot:
            plt.xlabel("z [%sm]"%scale)
            plt.ylabel("Amplitude [a.u.]")
           
        else:
            plt.axis('off')
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)

    if which in ["both", "phase"]:
        plt.sca(ax2)
        # ax2.set_aspect(1)
        pcm = plt.plot(screen.z/scale_factor, screen.phase[0,0,:])
    
        if not only_plot:
            plt.xlabel("z [%sm]"%scale)
            plt.ylabel("Phase [rad]")
                # plt.title("Phase")
        else:
            plt.axis('off')
            ax2.get_xaxis().set_visible(False)
            ax2.get_yaxis().set_visible(False) 
            







################## consider removing below
def plot_field_legacy(aperture, which='both', scale=None, colorbar=True, only_plot=False, filename=None, **kwargs):
    """
    Plots the given aperture

    Args:
        :aperture:  Aperture object of the mask
        :which:     default "both", "amplitude" or "phase"
        :colorbar:  True/False flag to plot colorbars
        :only_plot: if True, only shows image without labels and axes
        :filename:  if provided, saves figure to filename
        ###available output file extensions are same as opencv https://docs.opencv.org/3.4/d4/da8/group__imgcodecs.html
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
        plt.savefig("temp.png", bbox_inches='tight', pad_inches = 0)
        plt.close(fig1)
        img = cv2.imread("temp.png")
        cv2.imwrite(filename, img)

