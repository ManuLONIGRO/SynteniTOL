#!/bin/bash

cd /home/mlonigro/pvc/sintenia/ampliados_perfiles

# Lista de proteins unicas
#proteins=($(ls | awk -F'_' '{print $1}' | sort -u | grep -v "^profile$"))
# Lista de proteínas únicas a partir de los archivos *_cog
proteins=($(ls *_cog 2>/dev/null | sed 's/_cog$//' | sort -u))
cd ..

for protein in "${proteins[@]}";do
	/home/mlonigro/pvc/scripts/hmmbuild_noRepeat.sh "/home/mlonigro/pvc/sintenia/ampliados_perfiles/${protein}_cog"
	echo "Procesando ${protein}"
	/home/mlonigro/pvc/scripts/hmmsearch_noRepeat.sh "/home/mlonigro/pvc/sintenia/ampliados_perfiles/profile_${protein}.hmm"

done
