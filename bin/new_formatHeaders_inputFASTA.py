#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Nov 12 09:23 2025

Take the inputFASTA of NEXTFLOW and format the headers using dataframe.csv

@author: mlonigro
"""

import argparse
import sys
import pandas as pd
import os
from Bio import SeqIO
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input_fasta", type=str, required=True, help="Input FASTA file path")
parser.add_argument("--output_fasta", type=str, required=True, help="Output formatted FASTA file path")
parser.add_argument("--genomes_dir", type=str, required=False, help="Path to the directory containing .gbff files")
args = parser.parse_args()

input_fasta = args.input_fasta
output_fasta = args.output_fasta
genomes_dir = args.genomes_dir

def normalize(value):
    return (value or "NA").replace(" ", "_").replace(";","_").replace(":","_").replace("=","_").replace("(","_").replace(")","_")

# Get the protein_id from the fasta file. Format is >protein_id:other info
protein_ids = []
protein_id_to_seq = {}
with open(input_fasta, "r") as infile:
    for line in infile:
        if line.startswith(">"):
            # protein_id = line.split(':')[0].strip().replace(">","")
            header = line[1:].strip()
            match = re.match(r"^([^:\s]+):\d+-\d+(?:\s|$)", header)
            protein_id = match.group(1) if match else header.split(None,1)[0]  # Fallback to first word if regex fails
            protein_ids.append(protein_id)
    # Map each input protein_id to its full sequence (for the sequence-based fallback)
    infile.seek(0)
    input_protein_seqs = set()
    for record in SeqIO.parse(infile, "fasta"):
        rec_header = record.description
        match = re.match(r"^([^:\s]+):\d+-\d+(?:\s|$)", rec_header)
        rec_id = match.group(1) if match else rec_header.split(None,1)[0]
        seq = str(record.seq)
        protein_id_to_seq[rec_id] = seq
        input_protein_seqs.add(seq)

# Read gbff looking for taxid, species_name, isolate for each protein_id in protein_ids list
protein_id_to_info = {}
# Index of source-feature metadata per gbff, keyed by normalized organism (fallback)
genome_source_info = {}
for gbff_file in os.listdir(genomes_dir):
    if gbff_file.endswith(".gbff"):
        gbff_path = os.path.join(genomes_dir, gbff_file)
        for record in SeqIO.parse(gbff_path, "genbank"):
            organism = "NA"
            isolate = "NA"
            strain = "NA"
            taxid = "NA"
            for feature in record.features:
                if feature.type == "source":
                    db_xrefs = feature.qualifiers.get("db_xref", [])
                    for xref in db_xrefs:
                        if xref.startswith("taxon:"):
                            taxid = xref.split(":")[1]
                            break
                    organism = normalize(feature.qualifiers.get("organism", ["NA"])[0])
                    strain = normalize(feature.qualifiers.get("strain", ["NA"])[0])
                    isolate = normalize(feature.qualifiers.get("isolate", ["NA"])[0])

                # Search protein_id in CDS
                if feature.type == "CDS" and "protein_id" in feature.qualifiers:
                    pid = feature.qualifiers["protein_id"][0]
                    if pid in protein_ids:
                        protein_id_to_info[pid] = {
                            "taxid": (taxid if taxid is not None else "NA"),
                            "organism": (organism if organism is not None else "NA"),
                            "strain": (strain if strain is not None else "NA"),
                            "isolate": (isolate if isolate is not None else "NA"),                            
                        }
                    translation = feature.qualifiers.get("translation", [None])[0]
                    if translation and translation.strip() in input_protein_seqs:
                        key = f"{taxid}|{organism}|{strain}|{isolate}"
                        if key not in genome_source_info:
                            genome_source_info[key] = {
                                "taxid": taxid,
                                "organism": organism,
                                "strain": strain,
                                "isolate": isolate,
                                "seq": set(),
                            }
                        genome_source_info[key]["seq"].add(translation.strip())

# create header to replace for protein_id, create n_repetition of the header to add |n_repetitions to the new_header
header_count = {}
with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            # protein_id = line.split(':')[0].strip().replace(">","")
            header = line[1:].strip()
            match = re.match(r"^([^:\s]+):\d+-\d+(?:\s|$)", header)
            protein_id = match.group(1) if match else header.split(None,1)[0]  # Fallback to first word if regex fails
            
            if protein_id in protein_id_to_info:
                info = protein_id_to_info[protein_id]
                taxid = info.get("taxid", "")
                species_name = info.get("organism", "NA").replace(" ", "_").replace(";","_").replace(":","_").replace("=","_").replace("(","_").replace(")","_")
                strain = info.get("strain", "NA").replace(" ", "_").replace(";","_").replace(":","_").replace("=","_").replace("(","_").replace(")","_")
                isolate = info.get("isolate", "NA").replace(" ", "_").replace(";","_").replace(":","_").replace("=","_").replace("(","_").replace(")","_")

                # Count repetitions
                header = f">{taxid}|{species_name}|{strain}|{isolate}"
                if header not in header_count:
                    header_count[header] = 0
                header_count[header] += 1
                n_repetitions = header_count[header]
                new_header = f"{header}|{n_repetitions}\n"
                outfile.write(new_header)
            else:
                # Fallback: protein_id was not found as a CDS protein_id in any gbff
                # (e.g. RefSeq WP_ vs GenBank GCA_ assemblies). Try to assign the
                # source-feature metadata by matching the input protein sequence to
                # a CDS translation, so the tree header matches the iTOL datasets.
                seq = protein_id_to_seq.get(protein_id)
                if seq:
                    species_match = re.search(r"\[([^\]]+)\]", header)
                    species = normalize(species_match.group(1)) if species_match else None
                    seq_candidates = [info for info in genome_source_info.values() if seq in info["seq"]]
                    if species:
                        # Require the source organism to match the species in the header
                        candidates = [info for info in seq_candidates if info["organism"] == species or info["organism"].startswith(species + "_")]
                    else:
                        # No species available: only trust an unambiguous sequence match
                        candidates = seq_candidates if len(seq_candidates) == 1 else []
                    if candidates:
                        info = candidates[0]
                        taxid = info.get("taxid", "NA")
                        species_name = info.get("organism", "NA")
                        strain = info.get("strain", "NA")
                        isolate = info.get("isolate", "NA")

                        header = f">{taxid}|{species_name}|{strain}|{isolate}"
                        if header not in header_count:
                            header_count[header] = 0
                        header_count[header] += 1
                        new_header = f"{header}|{header_count[header]}\n"
                        outfile.write(new_header)
                    else:
                        print(f"WARNING: no reliable genome match for {protein_id}; keeping original header (header may not match iTOL datasets)", file=sys.stderr)
                        outfile.write(line)
                else:
                    print(f"WARNING: no sequence available for {protein_id}; keeping original header (header may not match iTOL datasets)", file=sys.stderr)
                    outfile.write(line)
        else:
            outfile.write(line)

