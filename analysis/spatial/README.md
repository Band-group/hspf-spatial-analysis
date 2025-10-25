# Snakemake workflow for `analysis/spatial`

This pipeline is organised into folders containing snakemake files that run various workflows. 

To ease the replication process a master file named **master.snakefile** calls snakefiles to run different modules or group of tasks as follows:
- [`rules/grid.snakefile`](rules/grid.snakefile): which generates global spatial polygons and grids
- [`rules/hbs.snakefile`](rules/hbs.snakefile): computes global HbS frequency surfaces using R-INLA.
- [`rules/hspf.snakefile`](rules/hspf.snakefile): implements the main regression models for *Pf* allele frequencies, as described in the paper
- [`rules/figures.snakefile`](rules/figures.snakefile): generates figure components for the manusript

The input HbS data for this pipeline can be found in [`input/cleanHbSdata.csv`](input/cleanHbSdata.csv) and the input *Pf* data can be found in [`input/hbs-pf-v8.sqlite`](input/hbs-pf-v8.sqlite).

There are also additional scripts which were used to prepare the HbS and Pf genotype data for the above analysis.

## Pipeline details
### Configuration

The master snakefile relies on a configuration file that specifies various aspects of what the pipeline will compute.  The config file used for the main analysis in the paper can be found in [`config/config-main.yaml`](config/config-main.yaml).  This includes choices on: the *Pf* dataset used, hyperparameters of the HbS model, the type and size of the grid cells used, a choice of covariates for HbS and *Pf* models, spatial criteria to select *Pf* data, the *Pf* loci of interest, and spatial domains of analysis. The file [`config/config-main.yaml`](config/config-full.yaml) includes a larger set of specifications required to replicate all results presented in the manuscript and supplementary material.

## Pipeline structure

Below we describe the major data and scripts for each snakemake workflow:
- **grid.snakefile**:
  - load world polygon: "geodata/naturalearthdata.Rdata"
  - generate spatial polygons: "code/create_aggregation_polygons.R"
- **extract_pf_data.snakefile**:
  - extract Pf data within spatial polygon: "input/scripts/extract_{dataset}_counts.R" (latest dataset: pf8)
- **hbs.snakefile**:
  - run HbS spatial model: "code/HbS_model_fit2.R"
  - plot fitted values from HbS model: "code/plot_HbS_fit.R"
  - load world polygon: "geodata/naturalearthdata.Rdata"
  - aggregate HbS data within spatial polygon: "code/aggregate_HbS_over_polygons.R"
  - aggregate covariate (raster) data within spatial polygon: "code/aggregate_raster_over_polygons.R"
  - compare model prediction of HbS model with Piel et al.: "code/plot_HbS_vs_piel_grid.R"
  - summarise fittend values from HbS model: "code/summarise_HbS_fits.R"
- **pf.snakefile**:
  - aggregate Pf values into polygons (longform format): "code/aggregate_pf_over_polygons_longform.R" )
  - aggregate ld values into polygons (longform format): "code/aggregate_pf_ld_over_polygons_longform.R" )
  - aggregate Pf (three ways ld) values into polygons (longform format): "code/aggregate_pf_3wayld_over_polygons_longform.R" )
- **hspf.snakefile**:
  - compile template model builder (TMB): "code/tmb/compile.R"
  - load world polygon: "geodata/naturalearthdata.Rdata"
  - aggregate covariate (raster) data within spatial polygon: "code/aggregate_raster_over_polygons.R"
  - run Pf model (BYM spatial regression): "code/BYM-tmb-longform.R"
  - plot fitted values from Pf momdel: "code/plot_hspf_fit.R"
  - map fitted values from Pf momdel: "code/plot_hspf_fit_grid.R"
  - summary of fitted values from Pf model: "code/summarise_hspf_fits.R"
  - put results of Pf model into supplementary table: "code/summary_hspf2excel.R"
- **figures.snakefile**:
  -  summary figure (manuscript): "code/figures/fig1.R",
  -  HbS-Pf association for various loci in Africa figure (manuscript): "code/figures/fig2_new.R",
  -  HbS-Pf time variation figure (manuscript): "code/figures/temporal_figure.R"
  -  linkage desequilibrium figure (manuscript): "code/figures/ld_figure.R"
  -  forest plot trend HbS-Pf (slope) figure (manuscript): "code/figures/forest_ggplot.R"
  -  data summary figure (manuscript): "code/data_summary.R" (optional)

## Pipeline installation

### Installing dependencies

The pipeline has several dependencies which must be installed before use.  Dependencies include:

- snakemake (v8 or above)
- R (we tested this using in version v4.4.0)
- The R-INLA package (tested using v22.08.24 and v24.12.11) and its dependencies.
- The TMB package (tested using v1.9.17)
- A standard set of packages for geographical modelling, including sf and stars.

Once R is installed, the function `install.prerequisites()` from [`code/functions.R`](code/functions.R) will attempt to identify and install all the needed R library prerequisites. Please make sure that the installation of all R packages was successful before running the pipeline.

### Data prerequisites

To replicate the analysis, you also need to create a `geodata/` folder (excluded from github) with the prerequisite large geographical data files in.  The needed prerequisite files are:

- HbS map by Piel et al.: `2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif` 
- spatial polygons (countries): `naturalearthdata.Rdata` 
- country borders (simplified): `ne_110m_admin_0_countries` ,`adm1`, and `adm2` folders containing shape files
- malaria prevalence (covariate) raster data: `2024_GBD2023_Global_PfPR_2000.tif`

A .zip file of the folder geodata (1.7 Gb) is available here:

[https://1drv.ms/f/s!At2csG1tWRgenchXzqpT8OK024wQLw?e=f16DTd](https://1drv.ms/u/c/1e18596d6db09cdd/EZk0TmliAO9MpTi_G0IP7d8B9y67py8Ht_uKT7sHGERo8Q?e=kyH35G)

Note that the code `geodata::elevation_global()` used in `functions.R`  will also download elevation data in the folder `geodata/`. The sources of the data are detailed in the paper.

Please contact the authors of the paper if you need help to access the data or install the prerequisite R packages.


