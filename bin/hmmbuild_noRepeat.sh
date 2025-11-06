#!/bin/bash

# Verificar si se ha proporcionado el nombre del archivo de entrada
if [ "$#" -ne 1 ]; then
    echo "Using: $0 file_cog"
    exit 1
fi

cd /home/mlonigro/pvc/sintenia/ampliados_perfiles/

# Asignar el nombre del archivo de entrada a una variable. El input es el archivo _cog
input_file="$1"

# Obtener el directorio del archivo de entrada
input_dir=$(dirname "$input_file")

# Obtener el nombre base del archivo de entrada (sin ruta y sin extensión)
base_name=$(basename "$input_file" _cog)

# Definir los nombres de los archivos de salida en el mismo directorio del archivo de entrada
cdhit_output="${input_dir}/${base_name}_hit"
mafft_output="${input_dir}/${base_name}_alin"
hmmbuild_output="${input_dir}/profile_${base_name}.hmm"

# Verificar si el modelo ya existe
if [ -f "$hmmbuild_output" ]; then
    echo "El modelo $hmmbuild_output ya existe. Saltando..."
    exit 0
fi

# Ejecutar cd-hit
cd-hit -i "$input_file" -o "$cdhit_output" -c 0.9

# Ejecutar mafft para alineamiento
mafft --auto --anysymbol "$cdhit_output" > "$mafft_output"

# Crear el perfil
hmmbuild "$hmmbuild_output" "$mafft_output"

echo "Proceso completado. $hmmbuild_output"
