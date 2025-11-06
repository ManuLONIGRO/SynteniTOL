#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 14:44:30 2025

Renombrar headers con el formato >taxid|species_name|gene_id
Crea ftr_blast_formated_species que contiene los nombres de las especies para itol 

@author: mlonigro
"""

# === 1. Crear índice protein_id → archivo.gbff ===
import os
from Bio import SeqIO

gbff_dir = "/home/mlonigro/pvc/sintenia/genomas"
index_file = "/home/mlonigro/pvc/sintenia/proteinid_gbff.json"  
protein_to_file = {}


print("Indexando archivos .gbff...")
for filename in os.listdir(gbff_dir):
    if filename.endswith(".gbff"):
        filepath = os.path.join(gbff_dir, filename)
        try:
            for record in SeqIO.parse(filepath, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS" and "protein_id" in feature.qualifiers:
                        pid = feature.qualifiers["protein_id"][0]
                        protein_to_file[pid] = filepath
        except Exception as e:
            print(f"⚠️ Error leyendo {filepath}: {e}")

import json
# Guardar índice para uso futuro (opcional)
with open(index_file, "w") as jf:
    json.dump(protein_to_file, jf)

print(f"Indexado completo: {len(protein_to_file)} proteinIDs indexados.\n")

#%%

import re
from ete3 import NCBITaxa
from Bio import SeqIO

import json

# Ruta al archivo que guardaste antes
index_file = "/home/mlonigro/pvc/sintenia/proteinid_gbff.json"

# Cargar el diccionario creado en la anterior celda
with open(index_file, "r") as jf:
    protein_to_file = json.load(jf)

# Inicializar NCBI Taxa
ncbi = NCBITaxa()

acc_in_fasta = []
dic_phyl = {}
genesID = []

#output_base_path = '/home/mlonigro/pvc/blast/fwd/'

subunit = "fwdc"

input_blast = f"/home/mlonigro/pvc/blast/09-10-2025/fwd/{subunit}_blast"
output_blast = f"/home/mlonigro/pvc/blast/09-10-2025/fwd/{subunit}_isolates"
genesID_file = f"/home/mlonigro/pvc/blast/09-10-2025/fwd/{subunit}_genesID.temporal"
protein_metadata = f"/home/mlonigro/pvc/blast/09-10-2025/fwd/{subunit}_prot_isolate.tsv.temporal"

# # PRUEBA
# input_blast = f"/home/mlonigro/pvc/clusterednr/fwd/{subunit}_blast"
# output_blast = f"/home/mlonigro/pvc/clusterednr/fwd/{subunit}_isolates"
# genesID_file = f"/home/mlonigro/pvc/clusterednr/fwd/{subunit}_genesID.temporal"
# protein_metadata = f"/home/mlonigro/pvc/clusterednr/fwd/{subunit}_prot_isolate.tsv.temporal"


# Leer archivo fasta del BLAST
#with open(input_blast) as infile, open("/home/mlonigro/pvc/blast/ftr/formated/fwd_blast_sinGenID2", "w") as outfile:
with open(input_blast, "r") as infile:
    for record in SeqIO.parse(infile, "fasta"):
        # Extraer protein ID (ej: WP_010962053.1)
        gene_id = record.id.split(":")[0]
        genesID.append(gene_id)
        # # Extraer especie entre corchetes: [Methylococcus capsulatus]
        # match = re.search(r"\[([^\]]+)\]", record.description)
        # if not match:
        #     continue
        # species_raw = match.group(1)
        # species_name = species_raw.replace(" ", "_")
        # species_clean = species_raw.strip()

        # # Buscar taxid con ete3
        # taxid = None
        # try:
        #     corrected_name, taxid = ncbi.get_fuzzy_name_translation(species_clean)
        # except:
        #     try:
        #         result = ncbi.get_name_translator([species_clean])
        #         if result:
        #             taxid = result[species_clean][0]
        #     except:
        #         pass

        # # Obtener phylum si se encontró taxid
        # phylum = "NOT_FOUND"
        # if taxid:
        #     try:
        #         lineage = ncbi.get_lineage(taxid)
        #         ranks = ncbi.get_rank(lineage)
        #         names = ncbi.get_taxid_translator(lineage)
        #         for tid in lineage:
        #             if ranks[tid] == "phylum": #phylum, kingdom, class
        #                 phylum = names[tid]
        #                 break
        #     except:
        #         pass
        # else:
        #     taxid = "NOT_FOUND"

        # # Guardar en diccionario
        # dic_phyl[f"{taxid}|{species_name}|{gene_id}"] = phylum # Para subunidades sindividuales
        # #dic_phyl[f"{taxid}|{species_name}"] = phylum # Para concatenaciones

        # # Reescribir header
        # #record.id = f"{taxid}|{species_name}|{gene_id}" #Puedo modificarlo para sacarge gene_id en caso de trabajar para hacer una concatenacion
        # record.id = f"{taxid}|{species_name}"
        # record.description = ""
        # #SeqIO.write(record, outfile, "fasta")

with open(genesID_file,"w") as genesIDfile:
    for i in genesID:
        genesIDfile.write(i + "\n")


#%
# === 2. Leer proteinIDs de interés ===
with open(genesID_file, "r") as f:
    target_ids = [line.strip() for line in f if line.strip()]

pid_problematics = []

# === 3. Buscar información por proteinID y guardar resultados ===
with open(protein_metadata, "w") as out:
    out.write("protein_id\tisolate\ttaxid\n")
    for pid in target_ids:
        isolate = "NA"
        taxid = "NA"

        filepath = protein_to_file.get(pid)
        if filepath:
            try:
                for record in SeqIO.parse(filepath, "genbank"):
                    found = False
                    for feature in record.features:
                        if feature.type == "source":
                            isolate = feature.qualifiers.get("isolate", ["NA"])[0].replace(" ","_")
                            for xref in feature.qualifiers.get("db_xref", []):
                                if xref.startswith("taxon:"):
                                    taxid = xref.split(":")[1]
                                    break
                            found = True
                            break
                    if found:
                        break
            except Exception as e:
                print(f"⚠️ Error procesando {filepath}: {e}")
        
        elif isolate == "NA" and taxid == "NA": #Guardar los pid que no fueron encontrados en ningun gbff
            pid_problematics.append(pid)

        out.write(f"{pid}\t{isolate}\t{taxid}\n")
        
print(f"✅ Proceso finalizado. Resultados en: {protein_metadata}")
#%
# pid_problematics_file = "/home/mlonigro/pvc/blast/fwd/pid_problematics"

# with open(pid_problematics_file, "w") as file:
#     for pid in pid_problematics:
#         file.write(pid + "\n")


#%
import re
from ete3 import NCBITaxa
from Bio import SeqIO

# Inicializar NCBI Taxa
ncbi = NCBITaxa()

acc_in_fasta = []
dic_phyl = {}
genesID = []
dic_isolate = {}
#output_base_path = '/home/mlonigro/pvc/blast/fwd/'

#protein_metadata = "/home/mlonigro/pvc/blast/fwd/protein_metadata.tsv"  # Se lo define en la anterior celda
with open (protein_metadata, "r") as protein_metadata_file:
    for line in protein_metadata_file:
        linea_spliteada = line.split("\t")
        protein_id = linea_spliteada[0]
        isolate_gbff = linea_spliteada[1]
        dic_isolate[protein_id] = isolate_gbff
        
# Leer archivo fasta del BLAST
with open(input_blast,"r") as infile, open(output_blast, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        # Extraer protein ID (ej: WP_010962053.1)
        gene_id = record.id.split(":")[0]
        genesID.append(gene_id)
        # Extraer especie entre corchetes: [Methylococcus capsulatus]
        match = re.search(r"\[([^\]]+)\]", record.description)
        if not match:
            continue
        species_raw = match.group(1)
        species_name = species_raw.replace(" ", "_")
        species_clean = species_raw.strip()

        # Buscar taxid con ete3
        taxid = None
        try:
            corrected_name, taxid = ncbi.get_fuzzy_name_translation(species_clean)
            species_name = corrected_name.replace(" ", "_")
        except:
            try:
                result = ncbi.get_name_translator([species_clean])
                if result:
                    taxid = result[species_clean][0]
                    name_dict = ncbi.get_taxid_translator([taxid])
                    if taxid in name_dict:
                        species_name = name_dict[taxid].replace(" ", "_")  # Usar nombre corregido
                
            except:
                pass

        # Obtener phylum si se encontró taxid
        phylum = "NOT_FOUND"
        if taxid:
            try:
                lineage = ncbi.get_lineage(taxid)
                ranks = ncbi.get_rank(lineage)
                names = ncbi.get_taxid_translator(lineage)
                for tid in lineage:
                    if ranks[tid] == "phylum": #phylum, kingdom, class
                        phylum = names[tid]
                        break
            except:
                pass
            
        else:
            taxid = "NOT_FOUND"

        # Relacionar con el diccionario de isolates
        isolate = dic_isolate[gene_id]

        # Guardar en diccionario
        #dic_phyl[f"{taxid}|{species_name}|{gene_id}"] = phylum # Para subunidades sindividuales
        dic_phyl[f"{taxid}|{species_name}|{isolate}"] = phylum # Para concatenaciones

        # Reescribir header
        #record.id = f"{taxid}|{species_name}|{gene_id}" #Puedo modificarlo para sacarge gene_id en caso de trabajar para hacer una concatenacion
        record.id = f"{taxid}|{species_name}|{isolate}"
        record.description = ""
        SeqIO.write(record, outfile, "fasta") # Escribir el archivo outfile

#with open("/home/mlonigro/pvc/blast/fwd/fwdc_genesID","w") as genesIDfile:
#    for i in genesID:
#        genesIDfile.write(i + "\n")
    