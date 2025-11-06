#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 10:41:17 2025

Create dataframe of synteny fwda - Wood Ljundahl pathway

UPDATED TO NEXTFLOW

Headers format:
taxid|species_name|isolate|n_repetitions

@author: mlonigro
"""


# Script para mapear IDs de proteínas encontradas por hmmsearch a sus locus_tag en archivos .gbff.
# Procesa múltiples archivos .result en paralelo para acelerar el análisis.
# Guarda los mapeos en 'id_locustag.txt' y los IDs no encontrados en 'ids_no_mapeadas.txt'.

import os
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Parametros
evalue_limit = 1e-5
formatted_evalue = f"{evalue_limit:.0e}" # Para expresarlo en notacion cientifica, para los nombres de los archivos
#coverage_threshold = 0.70
log_file = f"/home/mlonigro/pvc/sintenia/file2synteny_{formatted_evalue}.log"
dic_longitudes_output = f"/home/mlonigro/pvc/sintenia/locustag_longitud_{formatted_evalue}"
dic_aligned_output = f"/home/mlonigro/pvc/sintenia/locustag_aligned_{formatted_evalue}"
dic_evalue_output = f"/home/mlonigro/pvc/sintenia/locustag_evalue_{formatted_evalue}"

# Configuración de entrada
input_dir = "/home/mlonigro/pvc/sintenia/genomas/"
result_dir = "/home/mlonigro/pvc/sintenia/resultados_hmm_ampl/"
output_file = f"/home/mlonigro/pvc/sintenia/id_locustag_{formatted_evalue}.txt"
no_mapeadas_file = "/home/mlonigro/pvc/sintenia/ids_no_mapeadas.txt"

def procesar_archivo(filename):
    try:
        result_file = os.path.join(result_dir, filename)
        organism_name = filename.split("_")[0]
        gbff_file = os.path.join(input_dir, f"{organism_name}.gbff")
        name_file_gbff = f"{organism_name}.gbff"

        if not os.path.exists(gbff_file):
            return None

        with open(result_file, "r") as f:
            lines = f.readlines()

        if len(lines) >= 18 and "No hits detected" in lines[15]:
            return None

        hit_presence = None
        if len(lines) >= 15:
            columns = lines[14].split()
            if len(columns) > 8:
                evalue = float(columns[0])
                if evalue < evalue_limit:
                    hit_id = columns[8].split("|")[-1]
                    #recover_hit_id = hit_id
                    hit_presence = True
                    
        if not hit_presence:
            return None

        # To take into account in the coverage the sum of all the "domains" found.
        if hit_presence:

            # Search domains
            for i, line in enumerate(lines):
                if "Domain annotation for each sequence (and alignments):" in line:
                    line_hits = i + 4
                    break
            else:
                raise ValueError("No se encontró la sección de anotación de dominios\n")

            max_domains_to_check = 10
            total_aligned_length = 0
    
            for offset in range(max_domains_to_check):
                idx = line_hits + offset
                if idx >= len(lines):
                    break                
                
                fields = lines[idx].strip().split()
                
                if not fields or not fields[0].isdigit():
                    break  # End of the domain list
                
                try:
                    ali_from = int(fields[9])
                    ali_to = int(fields[10])
                    aligned_length = ali_to - ali_from + 1
                    total_aligned_length += aligned_length
                except Exception as e:
                    print("⚠️ Error procesando línea:")
                    print(lines[idx].strip())
                    print(e)

        if not hit_presence:
            return None
        if hit_presence:
            hit_id_clean = hit_id.split(".")[0]
            records = list(SeqIO.parse(gbff_file, "genbank"))
            for record in records:
                for feature in record.features:
                    if feature.type == "CDS":
                        qualifiers = feature.qualifiers
                        protein_id = qualifiers.get("protein_id", [""])[0]
                        locus_tag = qualifiers.get("locus_tag", [""])[0]
                        all_ids = [protein_id.split(".")[0]]# + [x.split(":")[-1] for x in db_xref]

                        if hit_id_clean in all_ids:# or locus_tag in all_ids: # if hit_id_clean in all_ids:
                            found_id = locus_tag if locus_tag else protein_id
                            sequence = feature.qualifiers["translation"][0]
                            length_seq = len(sequence)
                            if length_seq < 50:
                                print(f"Secuencia corta {locus_tag}: {sequence}")
                            return ("mapeado", hit_id_clean, found_id, length_seq, total_aligned_length,evalue)

            print(f"No se encontró {hit_id_clean} en {name_file_gbff}\n")
            return ("no_mapeado", hit_id_clean)
    except Exception as e:
        print(f"❌ Error procesando {filename}: {e}\n")
        return None 

if __name__ == "__main__":
    result_files = [f for f in os.listdir(result_dir) if f.endswith(".result")]
    
    n_cpus = multiprocessing.cpu_count()
    max_workers = max(1, int(n_cpus * 0.9))
    print(f"Usando {max_workers}/{n_cpus} núcleos disponibles")

    id_mapping = {}
    no_mapeadas = []
    dic_longitudes = {} #locus_tag: length_seq    
    dic_aligned = {} #locus_tag: numero de aa alineados.
    dic_evalue = {} # locus_tag: evalue

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(procesar_archivo, f): f for f in result_files}
        for future in as_completed(futures):
            resultado = future.result()
            if resultado:
                if resultado[0] == "mapeado":
                    _, seq_id, final_id,length_seq, total_aligned_length, evalue= resultado
                    id_mapping[seq_id] = final_id
                    dic_longitudes[final_id] = length_seq
                    dic_aligned[final_id] = total_aligned_length
                    dic_evalue[final_id] = evalue
                    
                elif resultado[0] == "no_mapeado":
                    _, hit_id_clean = resultado
                    no_mapeadas.append(hit_id_clean)

    # Guardar los mapeos
    with open(output_file, "w") as f:
        for seq_id, final_id in id_mapping.items(): # protein_id -> locus_tag
            f.write(f"{seq_id} -> {final_id}\n")

    # Guardar los no mapeados
    with open(no_mapeadas_file, "w") as f:
        for id_no_mapeado in no_mapeadas:
            f.write(f"{id_no_mapeado}\n")

    # Guardar el diccionario locus_tag: longitud de la secuencia
    with open(dic_longitudes_output,"w") as f:
        for final_id, length_seq in dic_longitudes.items():
            f.write(f"{final_id} -> {length_seq}\n")

    # Guardar el diccionario locus_tag: longitud de la secuencia alineada
    with open(dic_aligned_output,"w") as f:
        for final_id, total_aligned_length in dic_aligned.items():
            f.write(f"{final_id} -> {total_aligned_length}\n")

    # Guardar el diccionario locus_tag: evalue
    with open(dic_evalue_output,"w") as f:
        for final_id, evalue in dic_evalue.items():
            f.write(f"{final_id} -> {evalue}\n")

    # Mostrar resumen
    print(f"\n✅ Mapeos completados: {len(id_mapping)}")
    print(f"❌ No mapeados: {len(no_mapeadas)}\n")

with open(log_file, "w") as file:
    file.write(f"""El analisis fue realizado con un evalue de {formatted_evalue}\n
               Los archivos generados con ese evalue se guardaron en:\n
               {output_file}\n
               {no_mapeadas_file}\n
               {dic_longitudes_output}\n
               """)

#%%

"""
Obtengo el diccionario id_gene.txt que relaciona el ID que aparece en el nombre del archivo .gbff con 
el nombre del gen, para poder usar nombres del gen en el grafico en lugar de locus_tag
"""
import os

# Directorios
result_dir = "/home/mlonigro/pvc/sintenia/resultados_hmm_ampl/"
output_file = "/home/mlonigro/pvc/sintenia/id_gene.txt"

# Diccionario para almacenar el mapeo {hit_id: gene_name}
hmm_hit_mapping = {}

# Procesar archivos .result
for filename in os.listdir(result_dir):
    if filename.endswith(".result"):
        result_file = os.path.join(result_dir, filename)
        
        # Obtener el nombre del gen desde el nombre del archivo
        parts = filename.split("_")
        if len(parts) < 2:
            print(f"Nombre de archivo inesperado: {filename}")
            continue
        #gene_name = parts[1].replace(".result", "")  # Obtener el nombre del gen
        gene_name = parts[-1].replace(".result", "") # Obtener el nombre del gen que esta en la ultima parte al hacer el split
        # Leer el archivo y extraer el mejor hit (línea 17)
        with open(result_file, "r") as f:
            lines = f.readlines()

        if len(lines) >= 18 and "No hits detected" in lines[15]:
            print(f"No hits found in {result_file}")
            continue

        # Extraer el ID del mejor hit
        hit_id = None
        if len(lines) >= 15:
            columns = lines[14].split()
            if len(columns) > 8:
                hit_id = columns[8].split("|")[-1]# Extraer el ID del hit
                hit_id = hit_id.split(".")[0]

        if hit_id:
            hmm_hit_mapping[hit_id] = gene_name
        else:
            print(f"Error extrayendo hit en {result_file}")

# Guardar el mapeo en un archivo
with open(output_file, "w") as f:
    for hit_id, gene_name in hmm_hit_mapping.items():
        f.write(f"{hit_id} -> {gene_name}\n")

print(f"Proceso completado. Resultados guardados en {output_file}")



#%%
#-------------------------------------------------------------------
"Crear locustag_gene para relacionar los locus tag con los nombre de los genes y nombrarlos de esa forma en el grafico de la sintenia"

# Leer archivo 1: protein_id -> locus_tag
archivo1 = f"/home/mlonigro/pvc/sintenia/id_locustag_{formatted_evalue}.txt"
dict_locus = {}

with open(archivo1, "r") as f1:
    for line in f1:
        protein_id, locus_tag = line.strip().split(" -> ")  
        dict_locus[protein_id] = locus_tag

# Leer archivo 2: protein_id -> gene_name
archivo2 = "/home/mlonigro/pvc/sintenia/id_gene.txt"
dict_gene = {}

with open(archivo2, "r") as f2:
    for line in f2:
        protein_id, gene_name = line.strip().split(" -> ")  
        dict_gene[protein_id] = gene_name

# Relacionar locus_tag -> gene_name
locus_to_gene = {}

for protein_id, locus_tag in dict_locus.items():
    if protein_id in dict_gene:
        gene_name = dict_gene[protein_id]
        locus_to_gene[locus_tag] = gene_name

# Si quieres guardar el resultado en un archivo:
with open(f"/home/mlonigro/pvc/sintenia/locustag_gene_{formatted_evalue}", "w") as output_file:
    for locus_tag, gene_name in locus_to_gene.items():
        output_file.write(f"{locus_tag}\t{gene_name}\n")

print(f"Proceso completo guardado en: {archivo1}")

#%%
print("Graficando KDE...")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import seaborn as sns
import os
import math

# Lista ordenada deseada
gene_list = [
    "fwda", "fwdb", "fwdc", "ftr", "mch", "mtd", "hmd", "bmtd", "mer", "fae", "fdh", 
    "fold", "mcra", "mcrb", "mcrg", "pmoa", "pmob", "pmoc", "mfna", "mfnb", 
    "mfnc", "mfnd", "mfne", "mfnf", "mptd", "mpte", "mptg", "WP010870124", 
    "AAB98807.1", "AAB98766.1", "purb", "folp", "dmrx", "mptn", "cdha", "cdhb", "cdhc",
    "cdhd", "cdhe", "acsa", "acse", "ferredoxin", "acka","pta1","pta2","pta3",
]

# Crear carpeta para guardar la imagen compuesta
output_dir = f"/home/mlonigro/pvc/sintenia/densidades_por_proteina_{formatted_evalue}"
os.makedirs(output_dir, exist_ok=True)

#Cargar diccionario_longitud

# Agrupar longitudes por grupo
group_lengths = {}
for locus, group in locus_to_gene.items():
    length = dic_longitudes.get(locus)
    if length is not None:
        group_lengths.setdefault(group, []).append(length)

# Filtrar grupos válidos
valid_groups = {k: v for k, v in group_lengths.items() if len(v) >= 2}

# Determinar cuántos grupos se van a graficar (de los de gene_list que estén en valid_groups)
ordered_groups = [g for g in gene_list if g in valid_groups]
num_plots = len(ordered_groups)

# Crear figura
cols = 3
rows = math.ceil(num_plots / cols)
fig, axes = plt.subplots(rows, cols, figsize=(cols * 6, rows * 4))
axes = axes.flatten()

# Plotear según el orden de gene_list
for i, group in enumerate(ordered_groups):
    lengths = valid_groups[group]
    ax = axes[i]
    sns.kdeplot(lengths, fill=True, ax=ax)
    ax.set_title(f"Densidad - {group}")
    ax.set_xlabel("Longitud (aa)")
    ax.set_ylabel("Densidad")

# Ocultar ejes no usados
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.show()

# Guardar como imagen
formatted_evalue = f"{evalue_limit:.0e}"
output_path = os.path.join(output_dir, f"todas_las_densidades_{formatted_evalue}.png")
fig.savefig(output_path)
print(f"Imagen guardada en: {output_path}")

#%%
# Celda anterior pero multiproccessing
"""
presence_binary_data

ISOLATES

Junta todo los locustag de cada uno de los genes de interes de cada uno de los contigs en el mismo organismo
Sirve para graficar presencia/ausencia sin tener en cuenta la posicion.
hmm -> locus_tags -> gbff -> contexto genomico
"""
print("Inicio binary en multiproccessing")
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

# Leer archivo locus_tag -> gene_name
archivo1 = f"/home/mlonigro/pvc/sintenia/id_locustag_{formatted_evalue}.txt"
archivo2 = "/home/mlonigro/pvc/sintenia/id_gene.txt"

dict_locus = {}
dict_gene = {}

with open(archivo1, "r") as f1:
    for line in f1:
        protein_id, locus_tag = line.strip().split(" -> ")
        dict_locus[protein_id] = locus_tag

with open(archivo2, "r") as f2:
    for line in f2:
        protein_id, gene_name = line.strip().split(" -> ")
        dict_gene[protein_id] = gene_name

# Relacionar locus_tag -> gene_name
locus_to_gene = {dict_locus[protein_id]: dict_gene[protein_id] for protein_id in dict_locus if protein_id in dict_gene}

import os
from Bio import SeqIO
from ete3 import NCBITaxa
import multiprocessing as mp

ncbi = NCBITaxa()

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

    header = f"{taxid}|{species_name}|{isolate}"
    return (header, genome_length, all_genes)

# ------------------------------
# Código principal con Pool
# ------------------------------
def run_parallel(input_dir, locus_to_gene, n_threads=10):
    files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".gbff")]

    with mp.Pool(n_threads) as pool:
        results = pool.starmap(process_gbff, [(f, locus_to_gene) for f in files])

    # Filtrar los None
    genomes_data_binary = [r for r in results if r is not None]
    return genomes_data_binary

# Ejemplo de uso
input_dir = "/home/mlonigro/pvc/sintenia/genomas"
genomes_data_binary = run_parallel(input_dir, locus_to_gene, n_threads=10)
print(f"Se procesaron {len(genomes_data_binary)} genomas")

with open(f"/home/mlonigro/pvc/sintenia/presence_binary_data_{formatted_evalue}", "w") as f:
    for i in genomes_data_binary:
        f.write(str(i) + "\n")

#%%
print("Creando dataset taxid_locus_domain_full.tsv")

import os
#import matplotlib.pyplot as plt
#from Bio import SeqIO
from ete3 import NCBITaxa

choose_evalue = 1e-5
formatted_evalue = f"{choose_evalue:.0e}" # Para expresarlo en notacion cientifica, para los nombres de los archivos

# Inicializar NCBI Taxa
ncbi = NCBITaxa()

# --- Cargar genomes_data_binary ---
with open(f"/home/mlonigro/pvc/sintenia/presence_binary_data_{formatted_evalue}") as f:
    genomes_data_binary = [eval(line.strip()) for line in f]

# --- Caché para taxid -> (especie, dominio) ---
taxid_to_info = {}

def get_domain_and_species_name(taxid):
    if taxid in taxid_to_info:
        return taxid_to_info[taxid]
    try:
        lineage = ncbi.get_lineage(int(taxid))
        ranks   = ncbi.get_rank(lineage)
        names   = ncbi.get_taxid_translator(lineage)

        # Dominio
        domain_taxid = next((tid for tid in lineage if ranks[tid] == "domain"), None)
        domain = names.get(domain_taxid, "Desconocido") if domain_taxid else "Desconocido"

        # Especie
        species_taxid = next((tid for tid in lineage if ranks[tid] == "species"), None)
        species = names.get(species_taxid, "Desconocido") if species_taxid else "Desconocido"
        species = species.replace(' ','_')

    except Exception:
        domain, species = "Desconocido", "Desconocido"

    # Nota: guardamos en el orden (especie, dominio) para que coincida con el encabezado
    taxid_to_info[taxid] = (species, domain)
    return species, domain

# --- Función genérica para cargar mappings ---
def load_mapping(filepath, sep=" -> "):
    return {k: v for line in open(filepath) if sep in line for k, v in [line.strip().split(sep)]}

# --- Cargar diccionarios ---
id_to_locus     = load_mapping(f"/home/mlonigro/pvc/sintenia/id_locustag_{formatted_evalue}.txt")
id_to_gene      = load_mapping("/home/mlonigro/pvc/sintenia/id_gene.txt")
locus_to_length = load_mapping(f"/home/mlonigro/pvc/sintenia/locustag_longitud_{formatted_evalue}")
id_to_evalue    = load_mapping(f"/home/mlonigro/pvc/sintenia/locustag_evalue_{formatted_evalue}")
id_to_aligned   = load_mapping(f"/home/mlonigro/pvc/sintenia/locustag_aligned_{formatted_evalue}")

# Invertir diccionario
locus_to_id = {v: k for k, v in id_to_locus.items()}

dataframe = "/home/mlonigro/pvc/sintenia/dataframe_wlp.tsv"

# --- Generar archivo final ---
with open(dataframe, "w") as f_out:
    f_out.write("taxid\tespecie\tlocus_tag\tdominio\tid_proteina\tgene\tevalue\tlongitud\talineado\tcoverage\n")

    for header, _, genes in genomes_data_binary:
        taxid = header.split("|")[0]
        especie, dominio = get_domain_and_species_name(taxid)

        for *_, locus_tag in genes:
            id_proteina = locus_to_id.get(locus_tag, "NA")
            gene        = id_to_gene.get(id_proteina, "NA")
            longitud    = int(locus_to_length.get(locus_tag, "0"))
            aligned     = int(id_to_aligned.get(locus_tag, "0"))
            evalue      = id_to_evalue.get(locus_tag, "NA")
            coverage    = round(aligned / longitud, 2) if longitud > 0 else "NA"

            f_out.write(
                f"{taxid}\t{especie}\t{locus_tag}\t{dominio}\t{id_proteina}\t{gene}\t{evalue}\t{longitud}\t{aligned}\t{coverage}\n"
            )

print(f"Dataframe creado en {dataframe}")