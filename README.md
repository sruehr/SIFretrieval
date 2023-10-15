# SIFretrieval
The files in this repository are associated with **Ruehr et al. 2023. CITATION, DOI: XXXX**. The code provides basic functions and an example to process data from the Headwall Hyperspec Imager presented in the associated article. Four hyperspectral data cubes collected diurnally over one day (04-18-2022) are available at Zenodo repository **XXX**. The data provided in the repository are subsets of full data cubes and include only relevant spectral bands in order to reduce file size.

# Functions
The `functions.py` file includes functions to retrieve relevant indicies from the datacube produced by the hyperspectral imager. Solar-induced fluoresence (SIF) is retrieved using the 3 Fraunhofer Line Depth (3FLD) method with default bandwidths identified via the sensitivity analysis discussed in the article. Functions to retrieve the normalized difference vegetation index (NDVI), near-infrared reflectance of vegetation (NIRvR) and SIF yield ($\PhiF_sif$) are also included. 

# Example
The `example.py` file processes data provided in the Zenodo repository. Functions are loaded from the `functions.py` file to retrieve SIF, NDVI, NIRvR and Fy. Outputs are saved as separate files for further analysis (code not provided). 
