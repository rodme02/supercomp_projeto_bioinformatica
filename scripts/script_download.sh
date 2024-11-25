#!/bin/bash

# Criar uma pasta para os arquivos
mkdir -p dados
cd dados

# Loop para baixar e descompactar os arquivos de cromossomos
for X in {1..22}
do
    # Definir o URL de download
    URL="ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/snp147Mask/chr${X}.subst.fa.gz"
    
    # Baixar o arquivo
    echo "Baixando o cromossomo ${X}..."
    wget -q "$URL" -O "chr${X}.subst.fa.gz"
    
    # Descompactar o arquivo
    echo "Descompactando o cromossomo ${X}..."
    gunzip "chr${X}.subst.fa.gz"
done

echo "Download e descompactação concluídos."
