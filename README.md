# HbS-Pfsa spatial analysis repository

This repository contains code used for the analyses presented in

Python A. et al, "Geographical variation drives adaptive equilibrium of the P. falciparum sickle-associated mutations", https://doi.org/10.1101/2025.08.31.672853

It consists of

- scripts written in R (we have tested this using R v4.4.2)
- pipelines written in snakemake (we ran this using snakemake v9.3.0)
- C++ code for the TMB package.

The pipeline has been tested on CentOS Linux and on Mac OS X (Monterey and Sequoia).

See the enclosed LICENSE.txt file for the license.

## Installation

The pipeline has several dependencies which must be installed in R before use.  These include

- The R-INLA package (tested using v22.08.24 and v24.12.11) and its dependencies.
- The TMB package (tested using v1.9.17)
- The standard set of packages for geographical modelling, including sf and stars.
