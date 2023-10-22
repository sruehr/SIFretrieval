#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 09:54:30 2022

@author: sophieruehr

These functions are to retrieve and plot SIF and NDVI from the Headwall imager
"""

# Find index values in array closest to wavelegnths of interest 
def find_nearest(array, value):
    import numpy as np
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# Make colorbar plot 
def colplot(A, vmax = 1, vmin = 0, title = '',
            rois = None):
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ax = plt.subplot()
    ax.title.set_text(title)
    im = ax.imshow(A, vmax = vmax, vmin = vmin)
    divider = make_axes_locatable(ax)
    if rois != None:
         plt.scatter(x = rois[0],
                     y = rois[1], c = 'r')
         plt.scatter(x = rois[2],
                     y = rois[3], c = 'r')
        
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

# Retrieve NDVI using red (670) and infrared bands
def get_ndvi(wavel, lib, red = 670, IR = 800, bands = 5):
    import numpy as np
    red_index = find_nearest(wavel, red)
    IR_index = find_nearest(wavel, IR)
    R_band = lib[:,:,red_index:red_index+bands]
    IR_band = lib[:,:,IR_index-bands:IR_index]

    NDVI_denom = (R_band + IR_band)
    NDVI_num = (IR_band - R_band) 
    NDVI_denom = NDVI_denom.astype("float")
    NDVI_num = NDVI_num.astype("float")
    NDVI_denom[NDVI_denom == 0] = np.nan
    NDVI = NDVI_num / NDVI_denom
    NDVI = NDVI[:,:,0]
    return(NDVI)

# Retrieve NIR
def get_nir(wavel, lib, IR = 778, bands = 3):
    import numpy as np
    IR_index = find_nearest(wavel, IR)
    IR_band = lib[:,:,IR_index-bands:IR_index]
    IR = np.mean(IR_band, axis = 2)
    return(IR)

# Retrive RED
def get_red(wavel, lib, R = 671, bands = 3):
    import numpy as np
    R_index = find_nearest(wavel, R)
    R_band = lib[:,:,R_index-bands:R_index+bands]
    R = np.mean(R_band, axis = 2)
    return(R)

# Retrieve fPAR
def get_fpar(wavel, lib, roi_panel, RED = 671, bands = 5, roi_dim = 1):
    import numpy as np
    red_index = find_nearest(wavel, RED)
    red_band = lib[:,:,red_index-bands:red_index+bands]
    red_sum = np.sum(red_band, axis = 2)
    panel_red = roi_panel[:,:,red_index-bands:red_index+bands]
    panel_sum = np.sum(panel_red, axis = 2)
    
    ax = None
    if roi_dim == 2: # If normalizing row by row (white panel spans FOV)
        ax = 1
    panel_mean = np.mean(panel_sum, axis = ax)
    
    if roi_dim ==1:
        fPAR = 1-  red_sum / panel_mean
    if roi_dim == 2:
        fPAR = 1-  red_sum / panel_mean[:, None]

    return(fPAR)

# Take mean over wavelengths of interest for each pixel
def subset_lambda(lambda1, diff1, wavel1, lib1):
    import numpy as np
    lambdaL = lambda1 - diff1
    lambdaR = lambda1 + diff1
    index_L = find_nearest(wavel1, lambdaL)
    index_R = find_nearest(wavel1, lambdaR)
    lambda_calc = list(range(index_L,index_R,1))
    output = lib1[:,:,lambda_calc]
    output = np.nanmean(output, axis = 2)
    return(output)

# Retrieve SIF using 3 FLD method
# Default center wavelengths (lambdaL, lambdaR, lambdaIN) and bandwidths (diff_out, diff_in)
# derived from sensitivity analysis
def get_sif(wavel, roi_panel, lib, roi_dim = 1,
            lambdaL = 757.9163, lambdaR = 770.8011, lambdaIN = 760.6,
           diff_out = 0.4, diff_in = 0.4,
           mean_max_flag = 'mean', percents = [95, 5]):
        
    import numpy as np

    LoutL = subset_lambda(lambdaL, diff_out, wavel, lib)
    LoutR = subset_lambda(lambdaR, diff_out, wavel, lib)
    Lin = subset_lambda(lambdaIN, diff_in, wavel, lib)
    
    ax = None
    if roi_dim == 2:  # If normalizing row by row (white panel spans FOV)
        ax = 1
    
    # 3 methods to retrieve white panel values: mean, max, and percentiles.
    if mean_max_flag == 'mean':
        EoutL = np.nanmean(subset_lambda(lambdaL, diff_out, wavel, roi_panel), axis = ax)
        EoutR = np.nanmean(subset_lambda(lambdaR, diff_out, wavel, roi_panel), axis = ax)
        Ein = np.nanmean(subset_lambda(lambdaIN, diff_in, wavel, roi_panel), axis = ax)
        
    if mean_max_flag == 'max':
        EoutL = np.nanmax(subset_lambda(lambdaL, diff_out, wavel, roi_panel), axis = ax)
        EoutR = np.nanmax(subset_lambda(lambdaR, diff_out, wavel, roi_panel), axis = ax)
        Ein = np.nanmin(subset_lambda(lambdaIN, diff_in, wavel, roi_panel), axis = ax)
        
    if mean_max_flag == 'percent':
        EoutL = np.nanpercentile(subset_lambda(lambdaL, diff_out, wavel, roi_panel), percents[0], axis = ax)
        EoutR = np.nanpercentile(subset_lambda(lambdaR, diff_out, wavel, roi_panel), percents[0], axis = ax)
        Ein = np.nanpercentile(subset_lambda(lambdaIN, diff_in, wavel, roi_panel), percents[1], axis = ax)
        
        
    w_21 = (lambdaR - lambdaIN)/(lambdaR - lambdaL)
    w_22 = (lambdaIN - lambdaL)/(lambdaR - lambdaL)
    denom = 1 - (Ein/(w_21*EoutL+w_22*EoutR))

    if roi_dim == 2:
        Es = (Ein/(w_21*EoutL+w_22*EoutR))[:, None]
        numer = Lin - Es * (w_21*LoutL + w_22*LoutR)
        SIF = numer/denom[:, None]

    if roi_dim == 1:
        Es = (Ein/(w_21*EoutL+w_22*EoutR))
        numer = Lin - Es * (w_21*LoutL + w_22*LoutR)
        SIF = numer/denom

    return(SIF)


# Get Fyield following Zeng et al. 2019
def get_phif(SIF, NDVI, NIR):
    Fyield = SIF / (NDVI * NIR)
    return(Fyield)

# Spatillay upscale rasters (bilinear interpolation)
def resample_data(dat_path, upscale_factor):
    import rasterio
    from rasterio.enums import Resampling
    dataset = rasterio.open(dat_path)
    with rasterio.open(dat_path) as dataset:
    
        # resample data to target shape
        data = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.height * upscale_factor),
                int(dataset.width * upscale_factor)
            ),
            resampling=Resampling.bilinear
        )
    data = data[0,:,:]
    return(data)

