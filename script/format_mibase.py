import os
from Bio import SeqIO
import pandas as pd
import re

mirna_dat = snakemake.input[0]
out_mibase_file = snakemake.output[0]

with open(out_mibase_file, "w") as fout:
    fout.write("spe\tpre_id\tpre_seq\tmirna_id\tmirna_seq\tstart\tend\n")
    for record in SeqIO.parse(mirna_dat, "embl"):
        pre_id = record.name
        spe = pre_id[:3]
        pre_seq = record.seq.replace("U", "T")
        for ft in record.features:
            if "product" in ft.qualifiers:
                mirna_id = ft.qualifiers["product"][0]
            else:
                continue
            start = ft.location.start
            end = ft.location.end
            mirna_seq = pre_seq[start:end]
            fout.write(f"{spe}\t{pre_id}\t{pre_seq}\t{mirna_id}\t{mirna_seq}\t{start+1}\t{end}\n")
