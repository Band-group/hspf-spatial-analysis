This contains code and data for the spatial analysis.

## Changelog

** Update 17/03/2024 **

Moved the model-fitting code into `code/HbS_model_fit.R`, and the model diagnostic code into `code/HbS_model_diagnostics.R`.  Also added a pf regression diagnostic plot into the latter.  A few code tweaks have been made to make this work.

## Prerequisites

### R prerequisites

The function `install.prerequisites()` from `code/functions.R` will attempt to identify and install all the needed R library prerequisites. Good luck!

### Data prerequisites

Before starting you must create a `geodata/` folder (excluded from github) with the prerequisite large geographical data files in.  This includes:

- MAP malaria atlas data: `2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif` (162 Mb)
- (2012 version: `201201_Global_Sickle_Haemoglobin_HbS_Allele_Frequency_2010.tif`)
- Andre's 'naturalearth' R data files: `naturalearthdata.Rdata` (146 Mb)
- The `ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp` simplified country shape files.
- `pop100m.tif` which contains a processed population mask.
- `geodata::elevation_global()` will also download data to this folder (used in `functions.R` ).

For covariates we will additionally need:

- soil moisture Copernicus (nc files) (2.5 Gb)
- 2020_walking_only_travel_time_to_healthcare.tif (449 Mb)
- gpw-v4-population-density_2000.tif (213 Mb)

optional:

- `020_GBD2019_Global_PfPR_2019.tif`

For reference, geo data is being assembled here:

https://1drv.ms/f/s!At2csG1tWRgenchXzqpT8OK024wQLw?e=f16DTd

# README INPUT

## Storage limitations

To run Datapreparation.R you would need to obtain the following datasets:

- (optional for mapping only): 020_GBD2019_Global_PfPR_2019.tif (7.5 Mb)

To run HbS_Popmasking.R you would need to obtain the following datasets (population for various countries):

- fbookpopag.tif (1 Mb)
- sdn_ppp_2020.tif (897 Mb)
- ssd_ppp_2020.tif (311 Mb)
- som_ppp_2020.tif (329 Mb)
- eth_ppp_2020.tif (595 Mb)

Please send an email to the authors of the paper to access the data.

