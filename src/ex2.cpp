#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>
#include <filesystem>

namespace fs = std::filesystem;

using namespace std;

// Função para converter uma sequência de DNA para RNA
void convert_to_rna(const string& dna, string& rna) {
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < dna.size(); ++i) {
        if (dna[i] == 'T') {
            rna[i] = 'U';
        } else {
            rna[i] = dna[i];
        }
    }
}

// Função para processar um arquivo de DNA e gerar um arquivo de RNA
void process_file(const string& input_file, const string& output_file) {
    ifstream file(input_file);
    if (!file.is_open()) {
        cerr << "Erro ao abrir o arquivo: " << input_file << endl;
        return;
    }

    // Criar um arquivo para salvar o RNA
    ofstream output(output_file);
    if (!output.is_open()) {
        cerr << "Erro ao criar o arquivo de saída: " << output_file << endl;
        return;
    }

    // Processar o DNA linha por linha e converter para RNA
    string line;
    while (getline(file, line)) {
        string rna(line.size(), ' ');
        convert_to_rna(line, rna);
        output << rna << endl;
    }

    file.close();
    output.close();
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Configuração padrão: número de arquivos para processar
    int num_files = 22;
    if (argc > 1) {
        num_files = stoi(argv[1]);
    }

    // Lista de arquivos de entrada e saída
    vector<string> input_files, output_files;
    for (int i = 1; i <= num_files; ++i) {
        input_files.push_back("dados/processados/processed_chr" + to_string(i) + ".subst.fa");
        output_files.push_back("dados/rna/rna_chr" + to_string(i) + ".subst.fa");
    }

    // Criar diretório para salvar os arquivos de RNA, se não existir
    if (world_rank == 0 && !fs::exists("dados/rna")) {
        fs::create_directories("dados/rna");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Dividir arquivos entre processos
    int files_per_process = input_files.size() / world_size;
    int start_idx = world_rank * files_per_process;
    int end_idx = (world_rank == world_size - 1) ? input_files.size() : start_idx + files_per_process;

    // Processar os arquivos atribuídos a este processo
    for (int i = start_idx; i < end_idx; ++i) {
        cout << "Processo " << world_rank << " processando arquivo: " << input_files[i] << endl;
        process_file(input_files[i], output_files[i]);
    }

    MPI_Finalize();
    return 0;
}
