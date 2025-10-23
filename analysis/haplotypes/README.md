# Population-genetics analysis of MalariaGEN Pf7

This folder contains scripts for a population-genetic analysis of African samples from the MalariaGEN Pf7 dataset, obtained from:

https://www.malariagen.net/resource/34/

and described in doi:10.12688/wellcomeopenres.18681.1.

These analysis were conducted as part of

Python A. et al, "Geographical variation drives adaptive equilibrium of the P. falciparum sickle-associated mutations", doi:10.1101/2025.08.31.672853.


## Analysis structure:

The analysis is written as a snakemake pipeline in two parts:

- `pipeline/balancing/master.snakefile`

This pipeline conducts a process of data pre-preparation, including filtering data, removing/phasing mixed genotype calls, and polarising to the ancestral/derived status, followed by an analysis of selection metrics.

- `pipeline/relate/master.snakefile`

This pipeline implements the RELATE method to estimate genealogical trees and branch lengths across the genome.

The `scripts/` folder contains downstream analysis of these results, including:

- `scripts/haplotype_figure.R` which is the main haplotype figure in the paper
- `scripts/allele_age.R` implements a GEVA-style allele age estimate.
- `scripts/selection_metric_analysis.r` post-analyses the selection metrics, as reported in the paper.


