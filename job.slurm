#!/bin/bash
#SBATCH --job-name=mpi               # Nome do job
#SBATCH --output=mpi_%j.txt          # Nome do arquivo de saída
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=4096
#SBATCH --time=02:00:00
#SBATCH --partition=espec

mpirun -np 4 ./contar_bases