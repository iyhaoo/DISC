import loompy
import numpy as np
import pandas as pd
import collections
import argparse
import shutil
parser = argparse.ArgumentParser()

parser.add_argument('--input', required=True, type=str, help="loom")
FLAGS = vars(parser.parse_args())

loom_path = FLAGS["input"]
this_dir = loom_path.rsplit("/", 1)[0]
this_prefix = loom_path.rsplit(".", 1)[0]
output_path = "{}_unique_rename.loom".format(this_prefix)
with loompy.connect(loom_path, "r") as ds:
    gene_name = ds.ra["Gene"]
    print(gene_name)
gene_name_counter = collections.Counter(gene_name)
duplicated_gene_names = np.array(list(gene_name_counter.keys()))[np.array(list(gene_name_counter.values())) > 1]
output_message = "Has {} duplicated gene names".format(len(duplicated_gene_names))
print(output_message)
if len(duplicated_gene_names) == 0:
    with open("{}/{}".format(this_dir, output_message.replace(" ", "_")), "w"):
        pass
else:
    shutil.copy(loom_path, output_path)
    for gene_x in duplicated_gene_names:
        duplicated_gene_mask = gene_name == gene_x
        gene_name[duplicated_gene_mask] = np.array(["{}_duplicate_{}".format(x, y) for x, y in zip(gene_name[duplicated_gene_mask], range(duplicated_gene_mask.sum()))])
    print(pd.Series(gene_name).is_unique)
    with loompy.connect(output_path) as ds:
        ds.ra.Gene = gene_name
    with loompy.connect(output_path, "r") as ds:
        print(output_path, ds.ra["Gene"][duplicated_gene_mask])






