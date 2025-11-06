#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 10:52:59 2025

Creates the files for itol from dataframe.tsv

@author: mlonigro
"""
import os
from Bio import SeqIO
from ete3 import NCBITaxa
import multiprocessing as mp
from collections import defaultdict

ncbi = NCBITaxa()


import os
import matplotlib.pyplot as plt
#from matplotlib.patches import Polygon
from Bio import SeqIO
from ete3 import NCBITaxa

choose_evalue = 1e-5
formatted_evalue = f"{choose_evalue:.0e}" # Para expresarlo en notacion cientifica, para los nombres de los archivos

# Inicializar NCBI Taxa
ncbi = NCBITaxa()

# Directorio de entrada
input_dir = "/home/mlonigro/pvc/sintenia/genomas"
dataframe = "/home/mlonigro/pvc/sintenia/dataframe_wlp.tsv"

print(f"Inicio binary en multiproccessing from {dataframe}")

#--------------------------------------------------------------------------------------------------
# Cargar .tsv para filtrarlo y generar_fromTSV
import pandas as pd

# Leer el archivo .tsv
df = pd.read_csv(dataframe, sep="\t")

# Filtrar solo los que tengan coverage > 0.70
df_filtrado = df[df["coverage"] > 0.70]

# Crear diccionario {locus_tag: gen}
locus_to_gene = dict(zip(df_filtrado["locus_tag"], df_filtrado["gene"]))
dic_taxid = dict(zip(df_filtrado["locus_tag"],df_filtrado["taxid"]))
#---------------------------------------------------------------------------------------------------

print(f"Creating presence_binary_data_fromTSV_{formatted_evalue}...")

def process_gbff(input_file, locus_to_gene):
    records = list(SeqIO.parse(input_file, "genbank"))
    if not records:
        return None  

    all_genes = []
    genome_length = 0
    taxid = None
    isolate = "NA"

    for record in records:
        for feature in record.features:
            # Obtener taxid del source
            if feature.type == "source":
                db_xrefs = feature.qualifiers.get("db_xref", [])
                isolate = feature.qualifiers.get("isolate", ["NA"])[0].replace(" ", "_")
                for xref in db_xrefs:
                    if xref.startswith("taxon:"):
                        taxid = xref.split(":")[1]
                        break

            # Obtener genes de interés
            if feature.type == "CDS":
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                if locus_tag in locus_to_gene:
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = feature.location.strand
                    gene_name = locus_to_gene[locus_tag]
                    all_genes.append((start, end, strand, gene_name, locus_tag))

        if any(g[4] in locus_to_gene for g in all_genes):
            genome_length += len(record.seq)

    if not all_genes or not taxid:
        return None

    # Obtener nombre de especie
    try:
        lineage = ncbi.get_lineage(int(taxid))
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        species_taxid = [tid for tid in lineage if ranks[tid] == "species"]
        if species_taxid:
            nombre_org = names[species_taxid[0]]
            species_name = nombre_org.replace(" ", "_").replace("(", "").replace(")", "")
        else:
            species_name = "Unknown"
    except:
        species_name = "Unknown"

    header = f"{species_taxid[0]}|{species_name}|{isolate}" # Puse species_taxid porque en el gbff obtengo los de la cepa. Ver si coinciden con lso del arbol, o usar siempre el que encuentro y no elegir el del organismo o cepa.
    return (header, genome_length, all_genes)

# ------------------------------
# Código principal con Pool
# ------------------------------
def run_parallel(input_dir, locus_to_gene, n_threads=4):
    files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".gbff")]

    with mp.Pool(n_threads) as pool:
        results = pool.starmap(process_gbff, [(f, locus_to_gene) for f in files])

    # Filtrar los None
    genomes_data_binary = [r for r in results if r is not None]

    # Asignar numeración de repeticiones a los headers
    header_count = defaultdict(int)
    genomes_data_numbered = []
    for header, genome_length, all_genes in genomes_data_binary:
        header_count[header] += 1
        header_with_rep = f"{header}|{header_count[header]}"
        genomes_data_numbered.append((header_with_rep, genome_length, all_genes))
    
    return genomes_data_numbered

# Ejemplo de uso
input_dir = "/home/mlonigro/pvc/sintenia/genomas"
genomes_data_binary = run_parallel(input_dir, locus_to_gene, n_threads=8)
print(f"Se procesaron {len(genomes_data_binary)} genomas")

with open(f"/home/mlonigro/pvc/sintenia/presence_binary_data_fromTSV_{formatted_evalue}", "w") as f:
    for i in genomes_data_binary:
        f.write(str(i) + "\n")

#%%


#import os
#from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
# import random
import time
from ete3 import NCBITaxa

# Seleccionar gen

gene_of_interest = "fwda"

print(f"Creating genomic_context_data_isolates_{gene_of_interest}_{formatted_evalue} en MULTIPROCESSING")

# Inicializar NCBI Taxa
ncbi = NCBITaxa()

# ⏱️ Iniciar el temporizador
start_time = time.time()

# Directorio de entrada
input_dir = "/home/mlonigro/pvc/sintenia/genomas"


# Creo la lista manualmente porque necesito un orden específico para el archivo binario
gene_list = ["fwda","fwdb","fwdc","ftr","mch","mtd","hmd", "bmtd", "mer","fae","fdh","fold","mcra","mrcb","mcrg","pmoa","pmob","pmoc","mfna","mfnb","mfnc","mfnd","mfne","mfnf","mptd","mpte","WP010870124","AAB98807.1","CAA9891390.1","AAB98766.1","purb","folp","dmrx","mptn", "cdha", "cdhb", "cdhc","cdhd", "cdhe", "acsa", "acse", "ferredoxin","acka","pta1","pta2","pta3"]

header_count = {}

def process_file(filename):
    """Procesa un archivo .gbff y devuelve los datos si es válido."""
    input_file = os.path.join(input_dir, filename)
    records = list(SeqIO.parse(input_file, "genbank"))
    if not records:
        return None  # Saltar archivos sin datos

    contigs_with_protein = []
    genes = []
    isolate = "no_data"
    taxid = "no_data"
    found_metadata = False # Bandera para salir del bucle una vez que encuentre isolate y taxid

    # El contador de los headers, cuando son iguales va sumando 1 para el |1 |2 |3...
    header_count = Counter()
    results = []

    for record in records:
        contig_genes = []
        contains_gene = False  
        for feature in record.features:
            if feature.type == "source" and not found_metadata:
                db_xrefs = feature.qualifiers.get("db_xref", [])
                isolate = feature.qualifiers.get("isolate",["NA"])[0].replace(" ","_")
                for xref in db_xrefs:
                    if xref.startswith("taxon:"):
                        taxid = xref.split(":")[1]
                        break
                found_metadata = True
                
            if feature.type == "CDS":
                 qualifiers = feature.qualifiers
                 locus_tag = qualifiers.get("locus_tag", [""])[0]

                 if locus_tag in locus_to_gene:
                     start = int(feature.location.start)
                     end = int(feature.location.end)
                     strand = feature.location.strand
                     gene_name = locus_to_gene[locus_tag]
                     contig_genes.append((start, end, strand, gene_name, locus_tag))

                     if gene_of_interest in gene_name:
                         contains_gene = True
                 else:
                     start = int(feature.location.start)
                     end = int(feature.location.end)
                     strand = feature.location.strand
                     gene_name = locus_tag
                     contig_genes.append((start, end, strand, gene_name, locus_tag))

        if contains_gene:
            genes = contig_genes
            contigs_with_protein.append((record, contig_genes, len(record.seq)))
            if len(contigs_with_protein) > 1:
                print(f"{filename} contiene {len(contigs_with_protein)} contigos con {gene_of_interest}")

            # Elijo el que tenga mayor cantidad de genes alrededor del GOI
            record, genes, _ = max(contigs_with_protein, key=lambda x: len(x[1]))

    if not contigs_with_protein:
        return None  

    # Usar Ete3 para obtener el nombre del organismo con el taxid y generar el formato taxid|species_name
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    ranks = ncbi.get_rank(lineage)
    species_taxid = [tid for tid in lineage if ranks[tid] == "species"] # Elijo especie para no usar los taxids de cepas
    if species_taxid:   
        nombre_org = names[species_taxid[0]]
        species_name = nombre_org.replace(" ","_").replace("(","").replace(")","")
    else:
        print("No se pudo determinar el nombre de especie.")

    header = f"{species_taxid[0]}|{species_name}|{isolate}"
    header_count[header] += 1

    header_with_rep = f"{header}|{header_count[header]}"
    results.append((header_with_rep, len(record.seq), genes))  # suponiendo `genes` ya está definido
    
    return results

# Proceso principal: ejecutás en paralelo y juntás los resultados
genomes_data_synteny = []

# Uso de ProcessPoolExecutor para procesamiento en paralelo
with ProcessPoolExecutor() as executor:
    futures = {executor.submit(process_file, filename): filename for filename in os.listdir(input_dir) if filename.endswith(".gbff")}

    for future in as_completed(futures):
        result = future.result()
        if result and result not in genomes_data_synteny: # Agregue result not in genomes para evitar repeticiones
            genomes_data_synteny.append(result)

# Guardar resultados
with open(f"/home/mlonigro/pvc/sintenia/genomic_context_data_isolates_{gene_of_interest}_{formatted_evalue}", "w") as f:
    for item in genomes_data_synteny:
        f.write(str(item[0]) + "\n") # Cada item es una lista de un solo elemento.

# ⏱️ Detener el temporizador
end_time = time.time()
elapsed_time = end_time - start_time

# Imprimir el tiempo total de ejecución
print(f"\n🚀 Tiempo total de ejecución: {elapsed_time:.2f} segundos")

#%%

print("Creating reflect synteny...")
import ast

def reflect_synteny(genes):
    # Encontrar fwda
    fwda = next((g for g in genes if gene_of_interest in g[3].lower()), None)
    if not fwda:
        return genes  # No hay fwda, no se modifica nada

    fwda_start, fwda_end, fwda_strand = fwda[0], fwda[1], fwda[2]
    fwda_center = (fwda_start + fwda_end) // 2

    if fwda_strand == 1:
        return genes  # Ya está en la orientación deseada

    # Reflejar los genes alrededor de fwda_center
    reflected_genes = []
    for g in genes:
        start, end, strand, name, tag = g
        new_start = 2 * fwda_center - end
        new_end = 2 * fwda_center - start
        new_strand = -strand
        reflected_genes.append((min(new_start, new_end), max(new_start, new_end), new_strand, name, tag))

    # Ordenar por posición de inicio
    reflected_genes.sort(key=lambda x: x[0])
    return reflected_genes

# Leer el archivo original y escribir el orientado
input_path = f"/home/mlonigro/pvc/sintenia/genomic_context_data_isolates_{gene_of_interest}_{formatted_evalue}"
output_path = f"/home/mlonigro/pvc/sintenia/genomic_context_data_{gene_of_interest}_{formatted_evalue}_oriented"

with open(input_path, "r") as infile, open(output_path, "w") as outfile:
    for line in infile:
        data = ast.literal_eval(line.strip())
        header = data[0]
        genome_length = data[1]
        genes = data[2]

        genes_oriented = reflect_synteny(genes)
        outfile.write(f"{(header, genome_length, genes_oriented)}\n")

#%%
print("Creating context gene...")

import ast

def context_fwda(genes):
    # Encontrar fwda
    fwda = next((g for g in genes if gene_of_interest in g[3].lower()), None)
    if not fwda:
        return genes  # No hay fwda, no se modifica nada

    fwda_start, fwda_end, fwda_strand = fwda[0], fwda[1], fwda[2]
    fwda_center = (fwda_start + fwda_end) // 2
    
    
    fwda_context_upstream = fwda_center - 60000
    fwda_context_downstream = fwda_center + 60000
    
    # Reflejar los genes alrededor de fwda_center
    context_genes = []
    for g in genes:
        start, end, strand, name, tag = g
        if fwda_context_upstream < start < fwda_context_downstream: 
            context_genes.append(g)

    # Ordenar por posición de inicio
    context_genes.sort(key=lambda x: x[0])
    return context_genes

# Leer el archivo original y escribir el orientado
input_path = f"/home/mlonigro/pvc/sintenia/genomic_context_data_{gene_of_interest}_{formatted_evalue}_oriented"
output_path = f"/home/mlonigro/pvc/sintenia/genomic_context_data_{gene_of_interest}_{formatted_evalue}_sector"

with open(input_path, "r") as infile, open(output_path, "w") as outfile:
    for line in infile:
        data = ast.literal_eval(line.strip())
        header = data[0]
        genome_length = data[1]
        genes = data[2]

        genes_oriented = context_fwda(genes)
        outfile.write(f"{(header, genome_length, genes_oriented)}\n")