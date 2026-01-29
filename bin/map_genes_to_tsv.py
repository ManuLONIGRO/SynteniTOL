#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nextflow-friendly script to build the WLP dataframe TSV.

Inputs:
  --genomes_dir    Path to directory containing .gbff genomes
  --mappings_file  Path to merged mappings file (seq_id -> locus_tag\tlen=..\taligned=..\tevalue=..)
  --results_dir    Directory containing *.result files (to infer gene names)
  --out_tsv        Output TSV path

It reconstructs the dictionaries used by the original workflow and writes a
dataframe with columns:
  taxid, specie, locus_tag, domain, protein_id, gene, evalue, length, aligned coverage

New insight: created locus_to_coverage, with coverage from mappings file, "acc" column.
"""

import argparse
import os
from typing import Dict, Tuple, List
from Bio import SeqIO
# from ete3 import NCBITaxa

#from bin.createFiles2Synteny import locus_to_gene


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes_dir", required=True)
    parser.add_argument("--mappings_file", required=True)
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--out_tsv", required=True)
    return parser.parse_args()


# def load_id_to_gene_from_results(results_dir: str) -> Dict[str, str]:
#     id_to_gene: Dict[str, str] = {}
#     for name in os.listdir(results_dir):
#         if not name.endswith(".result"):
#             continue
#         # Gene name from filename, last token before .result
#         parts = name.split("_")
#         if len(parts) < 2:
#             continue
#         gene_name = parts[-1].replace(".result", "")

#         # Extract best-hit protein id from line 15 (index 14)
#         result_path = os.path.join(results_dir, name)
#         try:
#             with open(result_path, "r") as fh:
#                 lines = fh.readlines()
#             if len(lines) >= 18 and "No hits detected" in lines[15]:
#                 continue
#             if len(lines) >= 15:
#                 cols = lines[14].split()
#                 if len(cols) > 8:
#                     hit_id = cols[8].split("|")[-1].split(".")[0]
#                     id_to_gene[hit_id] = gene_name
#         except Exception:
#             continue
#     return id_to_gene


def load_mappings(mappings_file: str) -> Tuple[Dict[str, str], Dict[str, int], Dict[str, int], Dict[str, str], Dict[str, str], Dict[str, float]]:
    # Returns: id->locus, locus->len, locus->aligned, locus->evalue
    id_to_locus: Dict[str, str] = {}
    locus_to_length: Dict[str, int] = {}
    locus_to_aligned: Dict[str, int] = {}
    locus_to_evalue: Dict[str, str] = {}
    locus_to_gene: Dict[str, str] = {}
    locus_to_coverage: Dict[str, float] = {}

    with open(mappings_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            # Format: <seq_id> -> <locus>\tlen=<len>\taligned=<n>\tevalue=<e>
            try:
                left, right = line.split(" -> ", 1)
                seq_id = left
                fields = right.split("\t")
                locus_tag = fields[0]
                id_to_locus[seq_id] = locus_tag
                # Parse metrics if present
                for field in fields[1:]:
                    if field.startswith("len="):
                        try:
                            locus_to_length[locus_tag] = int(field.split("=", 1)[1])
                        except Exception:
                            pass
                    elif field.startswith("aligned="):
                        try:
                            locus_to_aligned[locus_tag] = int(field.split("=", 1)[1])
                        except Exception:
                            pass
                    elif field.startswith("evalue="):
                        try:
                            locus_to_evalue[locus_tag] = field.split("=", 1)[1]
                        except Exception:
                            pass
                    elif field.startswith("cog="):
                        try:
                            locus_to_gene[locus_tag] = field.split("=", 1)[1]
                        except Exception:
                            pass
                    elif field.startswith("hmmsearch_cover="):
                        try:
                            locus_to_coverage[locus_tag] = float(field.split("=", 1)[1])
                        except Exception:
                            pass
                    
            except ValueError:
                # line not conforming; skip
                continue

    return id_to_locus, locus_to_length, locus_to_aligned, locus_to_evalue, locus_to_gene, locus_to_coverage


# def get_tax_info_resolve():
#     ncbi = NCBITaxa()
#     cache: Dict[str, Tuple[str, str]] = {}

#     def resolve(taxid: str) -> Tuple[str, str]:
#         if taxid in cache:
#             return cache[taxid]
#         try:
#             lineage = ncbi.get_lineage(int(taxid))
#             ranks = ncbi.get_rank(lineage)
#             names = ncbi.get_taxid_translator(lineage)
#             domain_taxid = next((tid for tid in lineage if ranks[tid] == "domain"), None)
#             species_taxid = next((tid for tid in lineage if ranks[tid] == "species"), None)
#             phylum_taxid = next((tid for tid in lineage if ranks[tid] == "phylum"), None)
#             class_taxid = next((tid for tid in lineage if ranks[tid] == "class"), None)
#             domain = (names.get(domain_taxid, "Unknown") if domain_taxid else "Unknown").replace(" ", "_")
#             phylum = (names.get(phylum_taxid, "Unknown") if phylum_taxid else "Unknown").replace(" ", "_")
#             class_name = (names.get(class_taxid, "Unknown") if class_taxid else "Unknown").replace(" ", "_")
#             species = names.get(species_taxid, "Unknown") if species_taxid else "Unknown"
#             species = species.replace(" ", "_")

#         except Exception:
#             species, domain, phylum, class_name = "Unknown", "Unknown", "Unknown", "Unknown" 
#         cache[taxid] = (species, domain, phylum, class_name)
#            return species, domain, phylum, class_name
#     return resolve

def iterate_genomes(genomes_dir: str):
    for name in os.listdir(genomes_dir):
        if not name.endswith(".gbff"):
            continue
        yield os.path.join(genomes_dir, name)


def build_dataframe(genomes_dir: str,
                    id_to_locus: Dict[str, str],
                    # id_to_gene: Dict[str, str],
                    locus_to_length: Dict[str, int],
                    locus_to_aligned: Dict[str, int],
                    locus_to_evalue: Dict[str, str],
                    locus_to_gene: Dict[str,str],
                    locus_to_coverage: Dict[str,float],
                    out_tsv: str,
                    ) -> None:
    # resolve = get_tax_info_resolve()
    locus_to_id = {v: k for k, v in id_to_locus.items()}

    with open(out_tsv, "w") as out:
        # out.write("taxid\tdomain\tphylum\tclass\tspecie\tisolate\tlocus_tag\tprotein_id\tgene\tevalue\tlength\taligned\tcoverage\n")
        out.write("taxid\tlocus_tag\tprotein_id\tgene\tevalue\tlength\taligned\tcoverage\n")

        for gbff_path in iterate_genomes(genomes_dir):
            records = list(SeqIO.parse(gbff_path, "genbank"))
            if not records:
                continue
            taxid = None
            # local_isolate = "NA"
            for record in records:
                for feature in record.features:
                    if feature.type == "source":
                        db_xrefs = feature.qualifiers.get("db_xref", [])
                        # local_isolate = feature.qualifiers.get("isolate", ["NA"])[0].replace(" ", "_")
                        for xref in db_xrefs:
                            if xref.startswith("taxon:"):
                                taxid = xref.split(":")[1]
                                break
                        if taxid:
                            break
                if taxid:
                    break
            if not taxid:
                continue
            
            # specie, domain, phylum, class_name = resolve(taxid)

            # Emit rows for features whose locus_tag is in our mapping
            for record in records:
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    if not locus_tag:
                        continue
                    if locus_tag not in locus_to_id:
                        continue
                    protein_id = locus_to_id.get(locus_tag, "NA")
                    # gene = id_to_gene.get(protein_id, "NA")
                    gene = locus_to_gene.get(locus_tag, "NA")
                    length = int(locus_to_length.get(locus_tag, 0))
                    aligned = int(locus_to_aligned.get(locus_tag, 0))
                    evalue = float(locus_to_evalue.get(locus_tag, "NA"))
                    #coverage = round(aligned / length, 2) if length > 0 else "NA"
                    coverage = locus_to_coverage.get(locus_tag, "NA")

                    out.write(
                        # f"{taxid}\t{domain}\t{phylum}\t{class_name}\t{specie}\t{local_isolate}\t{locus_tag}\t{protein_id}\t{gene}\t{evalue}\t{length}\t{aligned}\t{coverage}\n"
                        f"{taxid}\t{locus_tag}\t{protein_id}\t{gene}\t{evalue}\t{length}\t{aligned}\t{coverage}\n"
                    )


def main():
    args = parse_args()
    id_to_locus, locus_to_length, locus_to_aligned, locus_to_evalue, locus_to_gene, locus_to_coverage = load_mappings(args.mappings_file)
    # id_to_gene = load_id_to_gene_from_results(args.results_dir)
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    build_dataframe(
        genomes_dir=args.genomes_dir,
        id_to_locus=id_to_locus,
        # id_to_gene=id_to_gene,
        locus_to_length=locus_to_length,
        locus_to_aligned=locus_to_aligned,
        locus_to_evalue=locus_to_evalue,
        locus_to_gene=locus_to_gene,
        locus_to_coverage=locus_to_coverage,
        out_tsv=args.out_tsv,
    )


if __name__ == "__main__":
    main()


