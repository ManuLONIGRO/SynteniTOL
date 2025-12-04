#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
To work on NEXTFLOW
Created on Tue Nov 25 10:51 2025
"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--protein_to_organism_map", required=True, help="Path to protein to organism mapping file")
parser.add_argument("--taxonomy_tsv", required=True, help="Path to taxonomy TSV file")
parser.add_argument("--out_pid_taxid_tsv", required=True, help="Path to output PID to TaxID TSV file")
args = parser.parse_args()

print("Loading tsv files ...")

df_prot = pd.read_csv(args.protein_to_organism_map, sep="\t")
df_tax = pd.read_csv(args.taxonomy_tsv, sep="\t")
print(f"Merging dataframes {df_prot} and {df_tax} ...")

df_out = df_prot.merge(df_tax, on="assembly", how= "left")
df_final = args.out_pid_taxid_tsv
df_out.to_csv(df_final, sep="\t", index=False)