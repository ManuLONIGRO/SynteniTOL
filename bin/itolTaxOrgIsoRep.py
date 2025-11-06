#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 09:59:04 2025

Crear itol files a partir de un concatenado o archivo ya formateado de alguna de estas formas:

1. taxid|species_name|geneID
2. taxid|species_name

ACTUALIZADO para NEXTFLOW. 22/10/2025

@author: mlonigro
"""

from Bio import SeqIO
from ete3 import NCBITaxa

ncbi = NCBITaxa()

# Diccionarios para almacenar phylum y dominio
dic_phyl = {}
dic_domain = {}
dic_clase = {}

import sys
from pathlib import Path

headers_file = Path(sys.argv[1])         # archivo generado por el primer proceso
output_dir = Path(sys.argv[2])           # carpeta de salida que define Nextflow
output_dir.mkdir(exist_ok=True, parents=True)

with open(headers_file) as infile:
    for line in infile:
        allparts = line.strip()
        parts = line.strip().split("|")
        # Variables por defecto
        taxid = "NOT_FOUND"
        species_name = "unknown"
        geneid = ""

        # Parseo del header
        taxid = parts[0]
        species_name = parts[1]
        isolate = parts[2]
        repetition = parts[3]

        # Inicializar valores
        phylum = "NOT_FOUND"
        domain = "NOT_FOUND"
        clase = "NOT_FOUND"
        
        # Si taxid es un número válido, consultar NCBI
        if taxid not in ("NOT_FOUND", "unknown"):
            try:
                taxid = int(taxid)
                lineage = ncbi.get_lineage(taxid)
                ranks = ncbi.get_rank(lineage)
                names = ncbi.get_taxid_translator(lineage)

                for tid in lineage:
                    if ranks.get(tid) == "phylum":
                        phylum = names.get(tid, "NOT_FOUND")
                        break
                for tid in lineage:
                    if ranks.get(tid) == "domain":
                        domain = names.get(tid, "NOT_FOUND")
                        break
                for tid in lineage:
                    if ranks.get(tid) == "class":
                        clase = names.get(tid, "NOT_FOUND")
                        break
            except Exception as e:
                print(f"⚠️ Error con taxid {taxid}: {e}")
                taxid = "NOT_FOUND"
                
        clave = allparts
        
        # Guardar en diccionarios
        dic_phyl[clave] = phylum
        dic_domain[clave] = domain
        dic_clase[clave] = clase
        
        

#-------------------------------------------------------------------------------------------------
# Valores unicos 
unique_values_phylum_conAsgard = list(set(dic_phyl.values())) #Contiene Asgard desglosado
unique_values_domain = list(set(dic_domain.values()))
unique_values_clase = list(set(dic_clase.values()))
#QUIERO reemplazar todos los phylums de Asgard por 1 solo nombre para no tener tantos en el arbol

# Lista de phylums que deben ser renombrados. IDEA: ESTOS PODRIAN SER CLASS
phylums_a_reemplazar = [
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
]

# Nuevo nombre
nuevo_phylum = "Promethearchaeati"

# Reemplazar los valores en el diccionario
for species, phylum in dic_phyl.items():
    if phylum in phylums_a_reemplazar:
        dic_phyl[species] = nuevo_phylum
    
    elif phylum == "NOT_FOUND":
        dic_phyl[species] = "Promethearchaeati"


unique_values_phylum = list(set(dic_phyl.values()))
#---------------------------------------------------------------------------------------
#phylum COLORES
import random

# Función para generar un color aleatorio en formato hexadecimal
def random_color():
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

# Generar un color para cada elemento en la lista
colors = [random_color() for _ in unique_values_phylum]
        
# Combinar elementos con sus colores correspondientes
elements_with_colors = list(zip(unique_values_phylum, colors))

list_colors_phylum = []

# Imprimir la lista de elementos con sus colores
for element, color in elements_with_colors:
    list_colors_phylum.append(color)
#--------------------------------------------------------------------------------------------
# domain COLORES

# Función para generar un color aleatorio en formato hexadecimal
def random_color():
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

# Generar un color para cada elemento en la lista
colors = [random_color() for _ in unique_values_domain]

# Combinar elementos con sus colores correspondientes
elements_with_colors = list(zip(unique_values_domain, colors))

list_colors_domain = []

# Imprimir la lista de elementos con sus colores
for element, color in elements_with_colors:
    list_colors_domain.append(color)
#--------------------------------------------------------------------------------------------
#clase COLORES

# Función para generar un color aleatorio en formato hexadecimal
def random_color():
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

# Generar un color para cada elemento en la lista
colors = [random_color() for _ in unique_values_clase]

# Combinar elementos con sus colores correspondientes
elements_with_colors = list(zip(unique_values_clase, colors))

list_colors_clase = []

# Imprimir la lista de elementos con sus colores
for element, color in elements_with_colors:
    list_colors_clase.append(color)
#--------------------------------------------------------------------------------------------



#Generar el itol_file
from tqdm import tqdm
from Bio import SeqIO

    # Especificar colores y etiquetas

phylum_labels = unique_values_phylum
phylum_colors = list_colors_phylum
phylum_color_map = dict(zip(phylum_labels, phylum_colors))

class_labels = unique_values_clase
class_colors = list_colors_clase
class_color_map = dict(zip(class_labels, class_colors))

domain_labels = unique_values_domain
domain_colors = list_colors_domain
domain_color_map = dict(zip(domain_labels, domain_colors)) 


    # Archivo para Phylum
with open(output_dir / "phylum.itol.txt", 'w') as phyl_file:
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
    for taxid, phylum in tqdm(dic_phyl.items()):
        taxid = taxid.replace("(","").replace(")","") # Los () generan problemas en el treefile, lo toma como un nuevo nodo.
        color = phylum_color_map.get(phylum, "#FFFFFF")  # Blanco si el phylum no está en el mapa
        phyl_file.write(f"{taxid},{color}\n")

    # Archivo para Clase
    with open(output_dir / "class.itol.txt", 'w') as class_file:
        class_file.write("DATASET_COLORSTRIP\n")
        class_file.write("SEPARATOR COMMA\n")
        class_file.write("DATASET_LABEL,class\n")
        class_file.write("COLOR,#ff0000\n")
        class_file.write("LEGEND_TITLE,Class\n")
        class_file.write("LEGEND_SHAPES," + ",".join(["1"] * len(class_labels)) + "\n")
        class_file.write("LEGEND_COLORS," + ",".join(class_colors) + "\n")
        class_file.write("LEGEND_LABELS," + ",".join(class_labels) + "\n")
        class_file.write("DATA\n")
        for taxid, class_ in tqdm(dic_clase.items()):
            taxid = taxid.replace("(","").replace(")","")
            color = class_color_map.get(class_, "#FFFFFF")  # Blanco si la clase no está en el mapa
            class_file.write(f"{taxid},{color}\n")

# Archivo para domain
with open(output_dir / "domain.itol.txt", 'w') as domain_file:
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
        taxid = taxid.replace("(","").replace(")","")
        color = domain_color_map.get(domain_, "#FFFFFF")  # Blanco si la clase no está en el mapa
        domain_file.write(f"{taxid},{color}\n")

# Generar archivos iTOL
#generate_itol_colorstrip_files(df, output_base_path)