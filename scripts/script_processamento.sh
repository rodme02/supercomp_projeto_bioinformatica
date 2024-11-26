#!/bin/bash

# Criar a pasta "processados" dentro de "dados"
mkdir -p dados/processados

# Loop por todos os arquivos na pasta "dados" com o padrão "chrX.subst.fa" onde X é de 1 a 22
for FILE in dados/chr*.subst.fa
do
    # Obter o nome do arquivo sem o caminho
    BASENAME=$(basename "$FILE")
    
    # Definir o nome do arquivo processado
    PROCESSED_FILE="dados/processados/processed_$BASENAME"
    
    # Remover a primeira linha, processar e dividir em linhas de 80 caracteres
    echo "Processando $BASENAME..."
    tail -n +2 "$FILE" | tr 'a-z' 'A-Z' | tr -d -c 'ATCG' | fold -w 80 > "$PROCESSED_FILE"
done

echo "Processamento concluído. Arquivos processados estão na pasta 'dados/processados/'."
