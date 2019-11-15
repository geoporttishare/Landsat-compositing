<img src="https://github.com/geoportti/Logos/blob/master/geoportti_logo_300px.png">

# Landsat compositing

This repository contains R scripts that can be used and edited for creating pixel-based image composites for Landsat TM/ETM+ data. Other Landsat sensors (MSS/OLI) have not yet been implemented. The scripts were adapted from a workflow described in "Van doninck J. and Tuomisto H. (2018) [A Landsat composite covering all Amazonia for applications in ecology and conservation]( https://doi.org/10.1002/rse2.77). Remote Sensing in Ecology and Conservation, 4(3), 197-210". 

## Contents

### 0. Wrapper script

The script **OGIIR_0_Wrapper_Example.R** provides a simple example on how the scripts can be used for downloading and processing Landsat TM/ETM+ data. Parameters set in this example script can be changed by the users to define their own spatial and temporal area of interest.

### 1. Query and download Landsat L1 data from USGS EarthExplorer
 
The script **OGIIR_1_browseEE.R** provides an interface to the [USGS/EROS Inventory Service Documentation](https://earthexplorer.usgs.gov/inventory/). Use of this services requires an [ERS](https://ers.cr.usgs.gov/login/) user profile authorization for Machine to Machine access, which must be requested separately ([custserv@usgs.gov](mailto:custserv@usgs.gov)).

### 2. Process L1 data to directionally normalized surface reflectance

The script **OGIIR_2_srTopoBRDF.R** uses Ledaps and other USGS-EROS scripts to convert Level-1 data to surface reflectance, to mask clouds and cloud cover, and to derive Landsat viewing angles. The script further provides options for normalization of topographic effects, if a digital elevation model is provided by the user, and for BRDF normalization.  

USGS-EROS ESPA product formatter ([https://github.com/USGS-EROS/espa-product-formatter](https://github.com/USGS-EROS/espa-product-formatter)) and ESPA Ledaps ([https://github.com/USGS-EROS/espa-surface-reflectance/tree/master/ledaps](https://github.com/USGS-EROS/espa-surface-reflectance/tree/master/ledaps)) must be installed to be able to run these scripts. 

### 3. Image projection and compositing

The scripts **OGIIR_3_composite.R** provides options for projecting a set of Landsat images to a common coordinate system and extent, and for performing pixel-based image compositing using a number of compositing criteria (maximum NDVI, band-wise median, multidimensional median).

## Citing
When used, the following citing should be mentioned: "We made use of geospatial data/instructions/computing resources provided by the Open Geospatial Information Infrastructure for Research (oGIIR, urn:nbn:fi:research-infras-2016072513) funded by the Academy of Finland.
