# Pseudoreplicate-Pipeline

# Running the Pipeline
The pipeline can either be run manually via command line or scheduled to run on a cluster via Slurm.

## Command Line
1. Open a terminal and navigate to the folder containing the workflow, e.g.
   
   `cd user/path_to_projects/Pseudoreplicate-Pipeline`

2. Activate conda environment with dependencies
3. Run the workflow via the command below, setting `n` as the number of processors to use

`snakemake --cores n`

## Running on a Cluster via Slurm
The pipeline can be scheduled to run on a cluster using the file `submit.sh`
