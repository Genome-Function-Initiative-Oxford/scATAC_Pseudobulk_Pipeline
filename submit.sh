#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=pseudobulk
#SBATCH --ntasks=8
#SBATCH --mem=100G
#SBATCH --time=01-00:00:00 # DAYS-HOURS:MINUTES:SECONDS
#SBATCH --output=slurm_debug/%j_%x.out
#SBATCH --error=slurm_debug/%j_%x.err

# Create debugging directory
mkdir slurm_debug -p

# Set source for conda environment (must be configured to match user directory)
source /path/to/baseenv/bin/activate pseudobulk

CONFIG_FILE=config/config.yaml

# Unlock workflow in case previous run failed
snakemake --configfile="$CONFIG_FILE" all --cores "$SLURM_NTASKS" --unlock
# Run workflow
snakemake --configfile="$CONFIG_FILE" all --cores "$SLURM_NTASKS" --rerun-incomplete --keep-going
