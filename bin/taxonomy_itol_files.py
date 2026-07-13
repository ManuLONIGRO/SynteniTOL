#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 09:59:04 2025

Copy pvc/scripts/itolTaxOrgIsoRep.py

UPDATED for NEXTFLOW 22/10/2026

@author: mlonigro
"""

import colorsys
from collections import defaultdict
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--taxonomy_tsv", required=True, help="Path to taxonomy.tsv file")
parser.add_argument("--df_complete_tsv", required=True, help="Path to df_complete.tsv file")
parser.add_argument("--fasta_file", required=True, help="Path to formatted_headers.fasta file")
parser.add_argument("--domain_itol_file", required=True, help="Output itol file for domain")
parser.add_argument("--phylum_itol_file", required=True, help="Output itol file for phylum")
parser.add_argument("--class_itol_file", required=True, help="Output itol file for class")

args = parser.parse_args()

# Load taxonomy.tsv
taxonomy_tsv = args.taxonomy_tsv

df_tax = pd.read_csv(taxonomy_tsv, sep="\t",names=["assembly","taxid","scientific_name","rank"], dtype={'assembly':str, 'taxid':str})
df_tax = df_tax[df_tax["rank"].isin(["domain","phylum","class"])]

# dictionary | assembly -> rank: scientific_name
assembly_to_taxonomy = defaultdict(dict)
for _, row in df_tax.iterrows():
    assembly_to_taxonomy[row["assembly"]][row["rank"]] = row["scientific_name"]

# Load df_complete.tsv
df_complete_tsv = args.df_complete_tsv
df_complete = pd.read_csv(df_complete_tsv, sep="\t",dtype={'local_taxid': str})

# Dictionary | taxid: assembly
taxid_to_assembly = dict(zip(df_complete["local_taxid"].astype(str), df_complete["assembly"]))

# Proccess fasta and build the dictionaries
from Bio import SeqIO

dic_domain = {}
dic_phylum = {}
dic_class = {}

fasta_file = args.fasta_file
for record in SeqIO.parse(fasta_file, "fasta"):
    header = record.id
    taxid = header.split("|")[0]
    if taxid not in taxid_to_assembly:
        continue
    assembly = taxid_to_assembly[taxid]

    if assembly not in assembly_to_taxonomy:
        continue

    tax = assembly_to_taxonomy[assembly]

    dic_domain[header] = tax.get("domain", "NA")
    dic_phylum[header] = tax.get("phylum", "NA")
    dic_class[header] = tax.get("class", "NA")

#-------------------------------------------------------------------------------------------------
# Valores unicos 
unique_values_phylum = list(set(dic_phylum.values())) #Contain desglosed groups of Asgard
unique_values_domain = list(set(dic_domain.values()))
unique_values_class = list(set(dic_class.values()))
# Replace all the phylums of Asgard by 1 single name to not have so many in the tree

# Phylums to be renamed
phylums_to_replace = [
    "Candidatus Heimdallarchaeota",
    "Candidatus Odinarchaeota",
    "Candidatus Thorarchaeota",
    "Candidatus Baldrarchaeota",
    "Candidatus Hermodarchaeota",
    "Candidatus Hodarchaeota",
    "Candidatus Helarchaeota",
    "Candidatus Borrarchaeota",
    "Candidatus Freyarchaeota", 
    "Promethearchaeota",
    "Candidatus Gerdarchaeota",
    "Candidatus Sigynarchaeota",
    "Candidatus Sifarchaeota",
    "Candidatus Wukongarchaeota",
    "Candidatus Kariarchaeota",
    "Candidatus Njordarchaeota",
    "Candidatus Idunnarchaeota",
    "Candidatus Friggarchaeota",
    "Candidatus Yidarchaeota",
]

# New name
new_phylum = "Promethearchaeati"

# Replace the values in the dictionary
for species, phylum in dic_phylum.items():
    if phylum in phylums_to_replace:
        dic_phylum[species] = new_phylum
    
    elif phylum == "NOT_FOUND":
        dic_phylum[species] = "Promethearchaeati"


unique_values_phylum_sinAsgard = list(set(dic_phylum.values()))


def evenly_spaced_colors(n):
    """High-contrast rainbow palette (same as bin/syntenyTaxOrg.py for COGs)."""
    colors = []
    for i in range(n):
        h = (i / max(1, n)) % 1.0
        s = 0.95
        l = 0.55
        r, g, b = colorsys.hls_to_rgb(h, l, s)
        colors.append(f"#{int(r * 255):02X}{int(g * 255):02X}{int(b * 255):02X}")
    return colors


def labels_and_colors(unique_values):
    labels = sorted(unique_values)
    colors = evenly_spaced_colors(len(labels))
    return labels, colors, dict(zip(labels, colors))


phylum_labels, phylum_colors, phylum_color_map = labels_and_colors(unique_values_phylum)
class_labels, class_colors, class_color_map = labels_and_colors(unique_values_class)
domain_labels, domain_colors, domain_color_map = labels_and_colors(unique_values_domain)

#Generate the itol_file
from tqdm import tqdm

# Output itol files
phylum_itol_file = args.phylum_itol_file
class_itol_file = args.class_itol_file
domain_itol_file = args.domain_itol_file

# File for Phylum
with open(phylum_itol_file, 'w') as phyl_file:
    phyl_file.write("DATASET_COLORSTRIP\n")
    phyl_file.write("SEPARATOR COMMA\n")
    phyl_file.write("COLOR_BRANCHES,1\n")
    phyl_file.write("DATASET_LABEL,phylum\n")
    phyl_file.write("COLOR,#ff0000\n")
    phyl_file.write("LEGEND_TITLE,Clades\n")
    phyl_file.write("LEGEND_SHAPES," + ",".join(["1"] * len(phylum_labels)) + "\n")
    phyl_file.write("LEGEND_COLORS," + ",".join(phylum_colors) + "\n")
    phyl_file.write("LEGEND_LABELS," + ",".join(phylum_labels) + "\n")
    phyl_file.write("DATA\n")
    for taxid, phylum in tqdm(dic_phylum.items()):
        taxid = taxid.replace("(","_").replace(")","_") # The "()" make troubles in the treefile, it take it like a new node
        color = phylum_color_map.get(phylum, "#FFFFFF")  # White if the phylum is not mapped
        phyl_file.write(f"{taxid},{color}\n")

# File for Class
with open(class_itol_file, 'w') as class_file:
    class_file.write("DATASET_COLORSTRIP\n")
    class_file.write("SEPARATOR COMMA\n")
    class_file.write("DATASET_LABEL,class\n")
    class_file.write("COLOR,#ff0000\n")
    class_file.write("LEGEND_TITLE,Class\n")
    class_file.write("LEGEND_SHAPES," + ",".join(["1"] * len(class_labels)) + "\n")
    class_file.write("LEGEND_COLORS," + ",".join(class_colors) + "\n")
    class_file.write("LEGEND_LABELS," + ",".join(class_labels) + "\n")
    class_file.write("DATA\n")
    for taxid, class_ in tqdm(dic_class.items()):
        taxid = taxid.replace("(","_").replace(")","_")
        color = class_color_map.get(class_, "#FFFFFF")  # White if the class is not mapped
        class_file.write(f"{taxid},{color}\n")

# File for domain
with open(domain_itol_file, 'w') as domain_file:
    domain_file.write("DATASET_COLORSTRIP\n")
    domain_file.write("SEPARATOR COMMA\n")
    domain_file.write("DATASET_LABEL,domain\n")
    domain_file.write("COLOR,#ff0000\n")
    domain_file.write("LEGEND_TITLE,domain\n")
    domain_file.write("LEGEND_SHAPES," + ",".join(["1"] * len(domain_labels)) + "\n")
    domain_file.write("LEGEND_COLORS," + ",".join(domain_colors) + "\n")
    domain_file.write("LEGEND_LABELS," + ",".join(domain_labels) + "\n")
    domain_file.write("DATA\n")
    for taxid,domain_ in tqdm(dic_domain.items()):
        taxid = taxid.replace("(","_").replace(")","_")
        color = domain_color_map.get(domain_, "#FFFFFF")  # White if the domain is not mapped
        domain_file.write(f"{taxid},{color}\n")