#!/bin/bash

# Verificar si el número adecuado de argumentos ha sido proporcionado
if [ $# -ne 1 ]; then
    #echo "Uso: $0 <archivo_perfil.hmm>"
    exit 1
fi

# Asignar el argumento proporcionado a una variable
profile_file=$1

# Verificar si el archivo de perfil existe
if [ ! -f "$profile_file" ]; then
    echo "Error: El archivo '$profile_file' no existe."
    exit 2
fi

# Crear directorio para resultados si no existe resultados_hmm_ampl
mkdir -p /home/mlonigro/pvc/sintenia/resultados_hmm_ampl/

#Extraer el nombre de la proteina a partir del archivo hmm
proteina=$(basename "$profile_file" .hmm | cut -d'_' -f2)

# Buscar todos los archivos .prot en el directorio proteomas
for proteoma in /home/mlonigro/pvc/sintenia/proteomas/*.prot; do
    # Extraer nombre base sin extensión
    nombre=$(basename "$proteoma" .prot)
    nombre=${nombre//_/}

    # Definir la ruta del archivo de resultados
    resultado="/home/mlonigro/pvc/sintenia/resultados_hmm_ampl/${nombre}_${proteina}.result"

    # Si el resultado ya existe, continuar al siguiente
    if [ -f "$resultado" ]; then
        echo "Ya existe resultado para $nombre, se omite."
        continue
    fi

    # Ejecutar hmmsearch
    #echo "Analizando $nombre..."
    hmmsearch "$profile_file" "$proteoma" > "$resultado"
done
