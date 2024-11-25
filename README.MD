# Projeto SuperComp - Bioinformática

### Rodrigo Paoliello de Medeiros

## Descrição

- `scripts/`: Scripts shell para baixar e processar os dados.
- `dados/`: Dados brutos e processados.
- `src/`: Códigos em C++ e executáveis.
- `slurms/`: Jobs para submissão no cluster.
- `outputs/`: Saídas dos jobs.

## Como rodar

1. Clone o repositório e entre na pasta do projeto:

```bash
git clone https://github.com/rodme02/supercomp_projeto_bioinformatica.git
cd supercomp_projeto_bioinformatica
```

2. Baixe os dados e faça o pré-processamento:

```bash
bash scripts/script_download.sh
bash scripts/script_processamento.sh
```

3. Compile cada exercício dessa forma:

```bash
mpic++ -fopenmp -o src/ex{x} src/ex{x}.cpp
```

4. Submeta os jobs de cada exercício dessa forma:

```bash
sbatch slurms/job_ex{x}.slurm
```

5. A saída de cada job estará na pasta `outputs/`.