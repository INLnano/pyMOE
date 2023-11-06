"""
plotting.py 
Module for plotting apertures ans save them as image files 


"""

import matplotlib.pyplot as plt
import cv2 

from pyMOE import Aperture
from pyMOE import Field
from pyMOE import ApertureField


def save_mask_plot(maskcir, xsize, ysize, filename):
    fig1 = plt.figure()
    figx = plt.imshow(maskcir, vmin=0, vmax=1,extent =[0,xsize,0,ysize], cmap=plt.get_cmap("Greys"))
    plt.axis('off')
    figx.axes.get_xaxis().set_visible(False)
    figx.axes.get_yaxis().set_visible(False)
    plt.savefig("temp.png", bbox_inches='tight', pad_inches = 0)
    plt.close(fig1)
    
    img = cv2.imread("temp.png")
    cv2.imwrite(filename, img)
    

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
        :field:  Aperture object of the mask
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

