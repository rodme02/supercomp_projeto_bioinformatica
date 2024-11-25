#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>
#include <filesystem>

namespace fs = std::filesystem;

using namespace std;

// Função para contar os códons de início (AUG) em uma sequência de RNA
int count_start_codons(const string& rna_sequence) {
    int count = 0;

    #pragma omp parallel for reduction(+:count)
    for (size_t i = 0; i <= rna_sequence.size() - 3; i += 3) {
        if (rna_sequence.substr(i, 3) == "AUG") {
            count++;
        }
    }

    return count;
}

// Função para processar um arquivo de RNA e contar os códons de início
int process_file(const string& file_path) {
    ifstream file(file_path);
    if (!file.is_open()) {
        cerr << "Erro ao abrir o arquivo: " << file_path << endl;
        return 0;
    }

    string line;
    int total_count = 0;

    while (getline(file, line)) {
        total_count += count_start_codons(line);
    }

    file.close();
    return total_count;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Número padrão de arquivos
    int num_files = 22;
    if (argc > 1) {
        num_files = stoi(argv[1]);
    }

    // Lista de arquivos RNA
    vector<string> rna_files;
    for (int i = 1; i <= num_files; ++i) {
        rna_files.push_back("dados/rna/rna_chr" + to_string(i) + ".subst.fa");
    }

    // Dividir arquivos entre processos
    int files_per_process = rna_files.size() / world_size;
    int start_idx = world_rank * files_per_process;
    int end_idx = (world_rank == world_size - 1) ? rna_files.size() : start_idx + files_per_process;

    // Contagem local de códons de início
    int local_count = 0;
    for (int i = start_idx; i < end_idx; ++i) {
        cout << "Processo " << world_rank << " processando arquivo: " << rna_files[i] << endl;
        local_count += process_file(rna_files[i]);
    }

    // Redução global para somar os resultados
    int global_count = 0;
    MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Apenas o processo de rank 0 exibe o resultado
    if (world_rank == 0) {
        cout << "Número total de proteínas inicializadas (AUG): " << global_count << endl;
    }

    MPI_Finalize();
    return 0;
}
