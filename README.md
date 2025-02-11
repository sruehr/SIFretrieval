# SIFretrieval
The files in this repository are associated with **Ruehr et al. 2023, "Quantifying seasonal and diurnal cycles of solar-induced fluorescence with a novel hyperspectral imager," a manuscript submitted to Geophysical Research Letters** ([https://doi.org/10.1029/2023GL107429](https://doi.org/10.1029/2023GL107429)). The code provides basic functions and an example to process data from the Headwall Hyperspec Imager presented in the associated article. Four hyperspectral data cubes collected diurnally over one day (04-18-2022) are available in this [repository ](https://zenodo.org/records/10246787). The data provided in the repository are subsets of full data cubes and include only relevant spectral bands in order to reduce file size.

## Functions for data processing
The `functions.py` file includes functions to retrieve relevant indicies from the datacube produced by the hyperspectral imager. Solar-induced fluoresence ($SIF$) is retrieved using the 3 Fraunhofer Line Depth (3FLD) method with default bandwidths identified via the sensitivity analysis discussed in the article. Functions to retrieve the normalized difference vegetation index ($NDVI$), near-infrared reflectance of vegetation ($NIRvR$) and a linear approximation of SIF yield ($\Phi F_{c}$) are also included. 

## Example of variable retrieval
The `example.py` file processes data provided in the Zenodo repository. Functions are loaded from the `functions.py` file to retrieve $SIF$, $NDVI$, $NIRvR$ and $\Phi F_{c}$. Outputs are saved as separate files for further analysis.
