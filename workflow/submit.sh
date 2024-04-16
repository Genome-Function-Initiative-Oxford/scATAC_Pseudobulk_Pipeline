#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=pseudo_snake
#SBATCH --ntasks=20
#SBATCH --mem=15G
#SBATCH --time=07-00:00:00 # DAYS-HOURS:MINUTES:SECONDS
#SBATCH --mail-user=tom.wilson@imm.ox.ac.uk
#SBATCH --output=slurm_debug/%j_%x.out
#SBATCH --error=slurm_debug/%j_%x.err

source /path/to/baseenv/bin/activate upstream

# select and edit the config file accordingly
snakemake --configfile=config/config.yaml all --cores 20 --unlock
snakemake --configfile=config/config.yaml all --cores 20 --rerun-incomplete #--keep-going