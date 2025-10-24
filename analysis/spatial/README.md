# Additional notes

### R prerequisites

In addition to the main readme file in the folder `hspf-spatial-analysis/analysis`, note that the function `install.prerequisites()` from `code/functions.R` will attempt to identify and install all the needed R library prerequisites. Please make sure that the installation of all R packages was successful.

### Data prerequisites

To replicate the analysis, you need to create a `geodata/` folder (excluded from github) with the prerequisite large geographical data files in.  This includes:

- HbS map: `2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif` 
- spatial polygons (countries): `naturalearthdata.Rdata` 
- Country borders (simplified): `ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp` 
- Population mask: `pop100m.tif` 
- Elevation data: `geodata::elevation_global()` will also download data to this folder (used in `functions.R` ).
- covariate raster data: `020_GBD2019_Global_PfPR_2019.tif`

For reference, geodata is being assembled here:

https://1drv.ms/f/s!At2csG1tWRgenchXzqpT8OK024wQLw?e=f16DTd

Please send an email to the authors of the paper if you need help to access the data.

