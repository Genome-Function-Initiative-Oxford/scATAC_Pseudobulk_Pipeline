#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=pseudo_snake
#SBATCH --ntasks=20
#SBATCH --mem=100G
#SBATCH --time=01-00:00:00 # DAYS-HOURS:MINUTES:SECONDS
#SBATCH --mail-user=tom.wilson@imm.ox.ac.uk
#SBATCH --output=slurm_debug/%j_%x.out
#SBATCH --error=slurm_debug/%j_%x.err

source /path/to/baseenv/bin/activate upstream

mkdir slurm_debug -p

# select and edit the config file accordingly
snakemake --configfile=config/config_2000.yaml all --cores 20 --unlock
snakemake --configfile=config/config_2000.yaml all --cores 20 --rerun-incomplete --keep-going