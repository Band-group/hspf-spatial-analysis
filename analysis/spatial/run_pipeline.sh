#!/bin/bash

# Update: I added --rerun-triggers mtime.
# This means only modification time is used to trigger pipeline reruns.
# In particular, code changes will not trigger re-runs!
args="$@"
snakemake \
--rerun-triggers mtime \
-s master.snakefile \
${args}

