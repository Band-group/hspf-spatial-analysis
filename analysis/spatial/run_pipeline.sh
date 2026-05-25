#!/bin/bash
args="$@"
snakemake \
-s master.snakefile \
${args}

