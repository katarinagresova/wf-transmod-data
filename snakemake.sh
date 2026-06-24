#!/bin/bash

conda activate snake

snakemake --use-conda -s workflow/Snakefile --cores 1