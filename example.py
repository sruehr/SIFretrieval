#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 17:35:06 2022

Author information:
Sophie Ruehr
Website: sruehr.github.io
Email: sophie.ruehr@berkeley.edu

Code used in Ruehr et al.; paper submitted to Geophysical Research Letters in November 2023

------------------------------------------------------------------------------
'Retrieve vegetation parameters from one hyperspectral acquisition'

This example code retrieves SIF, NDVI, NIR, fESC, and PhiF_c.
Retrievals are saved as .tif files
Functions for this code are stored in the `functions.py` file.
------------------------------------------------------------------------------
"""

#%% Input working directories -- change these to match your file system
# Paths to code, data, and save directories
    # Code (where function.py is stored)
code_wd = '../code'
    # Data
data_wd = '../data'
    # Save data (if desired, uncomment below)
save_wd = '../output'

# %% Load necessary libraries and functions 
import numpy as np
import glob
import spectral.io.envi as envi
import os
from tifffile import imsave
import scipy.ndimage
import matplotlib.pyplot as plt


# Import functions
os.chdir(code_wd)
from functions import colplot
from functions import get_sif
from functions import get_phif
from functions import get_ndvi
from functions import get_nir
from functions import get_red
from functions import find_nearest

#%% Define parameters
# Initial inputs
sif_function = 'percent' # Detect minimum and maximum values within 3FLD bands by 5 and 95 percentiles to reduce noise
ndvi_lim_percent = 25 # Quantile for NDVI mask below which vegetation is masked
upscale_factor = 1/4 # Bilinear resampling to larger spatial resolution 
file = '2022-04-18_11:20.hdr' # Header file for time of of interest; change if you want to process another time

# %% Open radiance datacube
# Set working directory to raw radiance data files
os.chdir(data_wd)
# Open header file (stores information/indices of hyperspectral datacube)
h = envi.read_envi_header(file)
# Retrieve information of wavelengths within datacube
wavel = [float(i) for i in h['wavelength']]
# Open radiance file 
lib = envi.open(file)

# %% Get panel radiance values for 3FLD SIF retrieval

# Plot first band of radiance data. Datacube is by default rotated 90˚ clockwise from horizontal
colplot(lib[:,:,0], vmax = 10)

# Identify ROI. Should be in form [x,y,x,y] (two points to define a rectangle)
roi = [750, 1500 ,800, 1600]
# Plot ROI on image to ensure accurate placement
colplot(lib[:,:,0], rois= roi, vmax=10)

# Extract white panel radiance 
roi_panel = lib[roi[1]:roi[3],roi[0]:roi[2],:]

# Plot spectra from one pixel
point = [100, 800]

# Plot spectra from 10 pixels. As this is a subset of the full data cube, there are missing values at some wavelengths
plt.plot(wavel, np.concatenate(np.concatenate(lib[0,point[1],:])))

# Zoom in to O2-band and its shoulders
lambdas = [757, 771]
inwavel = [find_nearest(wavel, lambdas[0]),find_nearest(wavel, lambdas[1])]
plt.plot(wavel[inwavel[0]:inwavel[1]], np.concatenate(np.concatenate(lib[0,point[1],inwavel[0]:inwavel[1]])))
                
# %% Retrieve SIF, NDVI, fPAR and PhiF
SIF = get_sif(wavel, roi_panel, lib, mean_max_flag = sif_function)
NDVI = get_ndvi(wavel, lib)
RED = get_red(wavel, lib) 
NIR = get_nir(wavel, lib)
             
# Rotate data 90 degrees (datacube is 90˚ rotated by default)
SIF = np.rot90(SIF) # 3FLD method
NDVI = np.rot90(NDVI)
NIR = np.rot90(NIR)
RED = np.rot90(RED)
            
# Calculate f_esc based on Zeng et al. 2019
PhiF = get_phif(SIF, NDVI, NIR)

# %% Resample data
# Spatially upscale data via bilinear interpolation to reduce noise and file size
NDVI = scipy.ndimage.zoom(NDVI, upscale_factor, order=1)
SIF = scipy.ndimage.zoom(SIF, upscale_factor, order=1)
NIR = scipy.ndimage.zoom(NIR, upscale_factor, order=1)
RED = scipy.ndimage.zoom(RED, upscale_factor, order=1)
PhiF =  scipy.ndimage.zoom(PhiF, upscale_factor, order=1)

# %% Plot output
colplot(SIF, vmax = 0.45, title = 'SIF')  
colplot(NDVI, title = 'NDVI') 
         
# %% Write data
os.chdir(save_wd)
imsave(('SIF.tif'), SIF) 
imsave(('NDVI.tif'), NDVI)
imsave(('NIR.tif'), NIR) 
imsave(('RED.tif'), RED) 
imsave(('PhiF.tif'), PhiF) 
    
