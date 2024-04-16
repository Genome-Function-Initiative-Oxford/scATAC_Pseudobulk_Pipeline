import pandas as pd
import os
import numpy as np

np.random.seed(1)

# Get parameters from config file
pseudo_size = snakemake.config["pseudo_size"]
min_replicates = snakemake.config["min_replicates"]
output_folder = snakemake.config["barcode_splitter_folder"]

# Input and output files
barcode_count_file = snakemake.input["barcode_count_file"]
cell_barcode_file = snakemake.output["cell_barcode_file"]

sample_name = snakemake.params["sample_name"]

# Create output folder if not yet created
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

try:
    # Open the cell barcode read counts as a DataFrame
    barcode_count_df = pd.read_csv(barcode_count_file, delim_whitespace = True, header = None, 
                                   names = ["Count", "CellBarcode"])
    print("Found", np.sum(barcode_count_df["Count"]), "reads across", len(barcode_count_df.index), 
          "cells for", sample_name)

    # Extract list of cell barcodes
    barcodes = np.array(list(barcode_count_df["CellBarcode"]))

except Exception as e:
    print("Could not open barcodes reads for", sample_name)
    raise(e)

if len(barcodes) == 0:
    raise Exception("No barcodes found for", sample_name)

# Calculate the number of pseudo-replicates that can be created each with n cells
n_pseudo_replicates = len(barcodes) // pseudo_size

if n_pseudo_replicates >= min_replicates:
    print("Splitting barcodes to create", n_pseudo_replicates, "pseudoreplicate(s) of size", 
          pseudo_size, "for", sample_name)
    # Randomly order barcodes
    np.random.shuffle(barcodes)

    # Create DataFrame with unique cell barcodes and assign them a pseudo-replicate number
    group_df = pd.DataFrame({"CellBarcode": barcodes[:(n_pseudo_replicates * pseudo_size)], 
                             "Group": np.repeat(np.array([sample_name + "_" + str(i + 1) for 
                                                          i in range(n_pseudo_replicates)]), pseudo_size)})
    # Save to file
    group_df.to_csv(cell_barcode_file, sep = '\t', index = False, header = False)

else:
    print("Not enough cells found to create", min_replicates, "or more pseudoreplicates for", sample_name)
    with open(cell_barcode_file, "w") as output_file:
        # Write message to file
        output_file.write("Found " + str(len(barcodes)) + " cell barcodes, which is not enough to create " +
                          str(min_replicates) + " of size " + str(pseudo_size) + " for " + sample_name)