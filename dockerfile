# Imagen base con Python y utilidades básicas
FROM ubuntu:22.04

# Evitar preguntas interactivas
ENV DEBIAN_FRONTEND=noninteractive

# Instalar dependencias básicas
RUN apt-get update && apt-get install -y \
    wget curl git build-essential openjdk-17-jre \
    hmmer cd-hit mafft fasttree \
    python3 python3-pip \
    && apt-get clean

# Instalar Entrez Direct (NCBI E-utilities)
RUN apt-get update && apt-get install -y wget curl unzip && \
    sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
ENV PATH="/root/edirect:${PATH}"

# Instalar Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
    && bash miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# Configurar canales, aceptar ToS y actualizar conda
RUN conda config --remove-key channels || true && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda config --set channel_priority strict && \
    conda update -y conda   

    
# Crear un entorno llamado syntenitol con Python 3.10 y librerías
RUN conda create -y -n syntenitol python=3.10

# Instalar herramientas adicionales desde bioconda
RUN conda install -y -n syntenitol bmge

# Instalar datasets de NCBI, ete3 biopython
RUN conda install -y -n syntenitol -c conda-forge ncbi-datasets-cli ete3 biopython pandas numpy

# Copiar scripts al contenedor
COPY bin /opt/syntenitol/bin

# Asegurarse que estan en el PATH
ENV PATH="/opt/syntenitol/bin:${PATH}"

# Definir el directorio de trabajo
WORKDIR /data