# scATAC Pseudobulk Pipeline
## About
This Snakemake workflow is designed to split one or more scATAC-seq BAM files into pseudobulk replicates, each containing n cells.

## Pipeline Steps
1. For each BAM file, count the number of unique occurences of each cell barcode
2. Assign cell barcodes a label corresponding to a pseudobulk replicate to create
3. Split the BAM file into pseudobulk BAM files using [Sinto](https://timoast.github.io/sinto/) to separate the cell barcodes
4. Generate indexes for the pseudobulk BAM files
5. Generate bigWigs for the pseudobulks
6. Call peaks for the pseudobulks
7. Create a metadata file summarising the number of pseudobulks created from each input BAM file

## Running the Pipeline
1. Set up the pseudobulk conda enviroment

   `conda env create --name pseudobulk --file=pseudobulk_environment.yml`
   
3. Edit config/config.yaml to set the data and parameters
4. Then the pipeline can either be run manually via command line or scheduled to run on a cluster via Slurm.

### Command Line
1. Open a terminal and navigate to the folder containing the workflow, e.g.
   
   `cd user/path_to_projects/scATAC_Pseudobulk_Pipeline`

2. Activate conda environment with dependencies
3. Run the workflow via the command below, setting `n` as the number of processors to use

`snakemake --cores n`

### Running on a Cluster via Slurm
The pipeline can be scheduled to run on a cluster using the file `submit.sh`
