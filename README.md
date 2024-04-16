# scATAC Pseudoreplicate Pipeline
## About
This Snakemake workflow is designed to split one or more scATAC-seq BAM files into pseudoreplicates, each containing n cells.

## Steps
1. For each BAM file, count the number of occurences for each unique cell barcode
2. Assign cell barcodes a label corresponding to a pseudoreplicate to create
3. Split the BAM file into pseudoreplicate BAM files using sinto to seperate the cell barcodes
4. Generate indexes for the pseudoreplicate BAM files
5. Generate bigWigs for the pseudoreplicates
6. Call peaks for the pseudoreplicates
7. Create a metadata file summarising the number of pseudoreplicates created from each input BAM file

## Running the Pipeline
The pipeline can either be run manually via command line or scheduled to run on a cluster via Slurm.

### Command Line
1. Open a terminal and navigate to the folder containing the workflow, e.g.
   
   `cd user/path_to_projects/Pseudoreplicate-Pipeline`

2. Activate conda environment with dependencies
3. Run the workflow via the command below, setting `n` as the number of processors to use

`snakemake --cores n`

### Running on a Cluster via Slurm
The pipeline can be scheduled to run on a cluster using the file `submit.sh`
