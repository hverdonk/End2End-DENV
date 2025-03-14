#!/bin/bash

snakemake \
-s snakefile-cfel-relax \
--printshellcmds \
--keep-going \
--rerun-incomplete \
--rerun-trigger mtime \
--cluster-cancel "qdel" \
--cluster "qsub -q lmem -V -l nodes=1:ppn={params.np} -l mem=256gb -o qsub/out/ -e qsub/err/ -l walltime=240:00:00" \
--jobs 10
