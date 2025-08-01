import os
import pandas as pd

# Get parameters from config file
pseudo_size = snakemake.config["pseudo_size"]

# Input and output files
finished_files = snakemake.input["finished_files"]
summary_csv = snakemake.output["summary_csv"]

file_prefix = snakemake.params["file_prefix"]

sample_names = []
n_pseudoreps = []
pseudorep_names = []

for file in finished_files:
    # Extract the sample name from the file path, e.g. "b-cell"
    sample_name = file.split(os.sep)[-1].split(".")[0].replace(file_prefix, "")

    with open(file) as txt:
        pseudoreps = [line.rstrip().split(os.sep)[-1].split(".")[0] for line in txt]

    sample_names.append(sample_name)
    if pseudoreps[0] != "":
        n_pseudoreps.append(len(pseudoreps))
        pseudorep_names.append(', '.join(pseudoreps))
    else:
        n_pseudoreps.append(0)
        pseudorep_names.append(None)

summary_df = pd.DataFrame({"Sample": sample_names, "nPseudoreps": n_pseudoreps, "PseudorepNames": pseudorep_names})
summary_df.to_csv(summary_csv, index = False, header = True)
