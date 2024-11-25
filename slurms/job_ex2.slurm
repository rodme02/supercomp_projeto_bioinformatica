#!/bin/bash
#SBATCH --job-name=ex2                       # Nome do job
#SBATCH --output=outputs/ex2_%j.txt          # Nome do arquivo de saída
#SBATCH --ntasks=4                           # Número de tarefas MPI  
#SBATCH --cpus-per-task=4                    # Número de threads OpenMP por tarefa
#SBATCH --mem=4096                                                     
#SBATCH --time=02:00:00
#SBATCH --partition=espec

mpirun -np 4 ./src/ex2 # Pode passar como argumento opcional a quantidade de arquivos a serem lidos (default = 22)