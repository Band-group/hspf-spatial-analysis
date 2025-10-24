# Additional notes

### R prerequisites

In addition to the main readme file in the folder `hspf-spatial-analysis/analysis`, note that the function `install.prerequisites()` from `code/functions.R` will attempt to identify and install all the needed R library prerequisites. Please make sure that the installation of all R packages was successful.

### Data prerequisites

To replicate the analysis, you need to create a `geodata/` folder (excluded from github) with the prerequisite large geographical data files in.  This includes:

- HbS map by Piel et al.: `2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif` 
- spatial polygons (countries): `naturalearthdata.Rdata` 
- country borders (simplified): `ne_110m_admin_0_countries` ,`adm1`, and `adm2` folders containing shape files
- malaria prevalence (covariate) raster data: `2024_GBD2023_Global_PfPR_2000.tif`

Note that the code `geodata::elevation_global()` used in `functions.R`  will download elevation data in the folder `geodata/`. The sources of the data are detailed in the paper.
For reference, a .zip file of the folder geodata (1.7 Gb) is available here:

[https://1drv.ms/f/s!At2csG1tWRgenchXzqpT8OK024wQLw?e=f16DTd](https://1drv.ms/u/c/1e18596d6db09cdd/EZk0TmliAO9MpTi_G0IP7d8B9y67py8Ht_uKT7sHGERo8Q?e=kyH35G)

Please contact the authors of the paper if you need help to access the data or install the prerequisite R packages.

