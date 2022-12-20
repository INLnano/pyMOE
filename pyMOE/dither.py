
"""
dither.py 
Module for dithering of masks/images 
"""

import cv2 
import numpy as np 

def floyd_steinberg(input_img , plot = False ): 
    """
    make dithered image 
    'input_img' = input 2D representation of img from cv2 
    optional: 
    plotting if True, shows the plot, deafults to False 
    """
    ###NOTE: considers the same pixel as in the image, possible improvement, change of pixel size 

    #input image as provided 
    img_gray_eq = input_img
    
    h,w = img_gray_eq.shape

    img_dither = np.zeros((h+1, w+1), dtype=np.float)
    
    threshold = 128

    for i in np.arange(0,h):
        for j in np.arange(0,w):
            img_dither[i, j] = img_gray_eq[i, j]
    
    for i in np.arange(0,h):
        for j in np.arange(0,w):
            opix = img_dither[i, j]
            if (img_dither[i, j] > threshold):
                vpix = 255
            else:
                vpix = 0

            img_dither[i, j] = vpix

            err = opix - vpix

            if j > 0:
                img_dither[i+1, j-1] = img_dither[i+1, j-1] + err * 3 / 16
            img_dither[i+1, j] = img_dither[i+1, j] + err * 5 / 16
            img_dither[i, j+1] = img_dither[i, j+1] + err * 7 / 16
            img_dither[i+1, j+1] = img_dither[i+1, j+1] + err * 1 / 16

    img_dither = img_dither.astype(np.uint8)
    img_dither = img_dither[0:h, 0:w]
    
    if plot == True: 
        import matplotlib.pyplot as plt 
        fig = plt.figure()
        plt.imshow(img_dither, vmin=0, vmax=255, cmap=plt.get_cmap("Greys"))

    #cv2.imwrite(output_filename, img_dither_inv)
    
    return img_gray_eq, img_dither
    
    
def dither_img(input_img, output_filename, plotting = False): 
    """
    make dithered image from input_img, save to output_img
    'inputcv2' = input 2D representation of img from cv2 
    'output_filename' = filename of image to be written 
    'pl' = plotting parameter, defaults to False 
    """

    img_gray0 = cv2.imread(input_img, cv2.IMREAD_GRAYSCALE)
    img_gray0 = 255 - img_gray0
    
    #possible IMPROVEMENT considering passing any dithering algorithm as argument to to dither_img 
    img_gray_eq, img_dither= floyd_steinberg(img_gray0, plot = plotting)

    cv2.imwrite(output_filename, img_dither)
    
    del img_gray_eq, img_dither
    
    
    
####missing other possible algorithms 



