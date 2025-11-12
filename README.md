# HbS-Pfsa spatial analysis repository

This repository contains code used for the analyses presented in

Python A. et al, "Geographical variation drives adaptive equilibrium of the P. falciparum sickle-associated mutations", https://doi.org/10.1101/2025.08.31.672853

It consists of a set of pipelines written in [snakemake](https://snakemake.readthedocs.io/en/stable/) which in turn use a set of [R](https://cran.r-project.org) scripts to conduct analysis.  There is also some C++ code used for spatial model fitting, written using [TMB](https://kaskr.github.io/adcomp/Introduction.html).

The pipelines have been run on CentOS Linux and on Mac OS X (Monterey and Sequoia).

In the [`theory/html/hspf-gpu`](./theory/html/hspf-gpu) folder there is also a javascript/WebGPU implementation of the parasite evolution model described in the above paper.

This code is licensed using the BOOST software license - see the enclosed [LICENSE.txt](./LICENSE.txt) file for details.

## Structure

Please see the README files in the following subfolders for more information, including installation and prerequisites:

- The [`analysis/spatial`](./analysis/spatial) folder contains the main pipeline for geospatial analyses used in the paper.
- The [`analysis/haplotypes`](./analysis/haplotypes) folder contains analysis of natural selection, haplotype and genealogical structure.
- The [`theory/html/hspf-gpu`](./theory/html/hspf-gpu) folder contains javascript and WebGPU code implementing an interactive simulation of parasite evolution over the HbS allele frequency surface.

