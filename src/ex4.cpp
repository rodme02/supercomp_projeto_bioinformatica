#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <mpi.h>
#include <omp.h>

using namespace std;

// Mapeamento de códons para aminoácidos (representados por números)
unordered_map<string, int> codon_to_aminoacid = {
    {"CCA", 1}, {"CCG", 1}, {"CCU", 1}, {"CCC", 1}, // Prolina
    {"UCU", 2}, {"UCA", 2}, {"UCG", 2}, {"UCC", 2}, // Serina
    {"CAG", 3}, {"CAA", 3},                         // Glutamina
    {"ACA", 4}, {"ACC", 4}, {"ACU", 4}, {"ACG", 4}, // Treonina
    {"AUG", 5},                                     // Metionina (Início)
    {"UGC", 6}, {"UGU", 6},                         // Cisteína
    {"GUG", 7}, {"GUA", 7}, {"GUC", 7}, {"GUU", 7}, // Valina
    {"UGA", -1}                                     // Códon STOP
};

// Função para traduzir uma sequência de RNA em aminoácidos até o códon de parada
vector<int> translate_rna(const string& rna_sequence) {
    vector<int> protein;
    #pragma omp parallel
    {
        vector<int> local_protein;
        #pragma omp for schedule(static)
        for (size_t i = 0; i <= rna_sequence.size() - 3; i += 3) {
            string codon = rna_sequence.substr(i, 3);
            if (codon_to_aminoacid.find(codon) != codon_to_aminoacid.end()) {
                int aminoacid = codon_to_aminoacid[codon];
                if (aminoacid == -1) break; // STOP codon found
                local_protein.push_back(aminoacid);
            }
        }
        #pragma omp critical
        protein.insert(protein.end(), local_protein.begin(), local_protein.end());
    }
    return protein;
}

// Função para processar um arquivo de RNA
vector<int> process_file(const string& file_path) {
    ifstream file(file_path);
    if (!file.is_open()) {
        cerr << "Erro ao abrir o arquivo: " << file_path << endl;
        return {};
    }

    vector<int> full_protein;
    string line;

    while (getline(file, line)) {
        vector<int> protein = translate_rna(line);
        full_protein.insert(full_protein.end(), protein.begin(), protein.end());
    }

    file.close();
    return full_protein;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Número padrão de arquivos RNA
    int num_files = 22;
    if (argc > 1) {
        num_files = stoi(argv[1]);
    }

    // Lista de arquivos de RNA
    vector<string> rna_files;
    for (int i = 1; i <= num_files; ++i) {
        rna_files.push_back("dados/rna/rna_chr" + to_string(i) + ".subst.fa");
    }

    // Dividir arquivos entre processos
    int files_per_process = rna_files.size() / world_size;
    int start_idx = world_rank * files_per_process;
    int end_idx = (world_rank == world_size - 1) ? rna_files.size() : start_idx + files_per_process;

    // Contagem local de proteínas traduzidas
    vector<int> local_protein;
    for (int i = start_idx; i < end_idx; ++i) {
        cout << "Processo " << world_rank << " processando arquivo: " << rna_files[i] << endl;
        vector<int> protein = process_file(rna_files[i]);
        local_protein.insert(local_protein.end(), protein.begin(), protein.end());
    }

    // Redução global para consolidar as proteínas traduzidas
    vector<int> global_protein;
    if (world_rank == 0) {
        global_protein.resize(local_protein.size() * world_size);
    }

    MPI_Gather(local_protein.data(), local_protein.size(), MPI_INT,
               global_protein.data(), local_protein.size(), MPI_INT, 0, MPI_COMM_WORLD);

    // Apenas o processo mestre exibe o resultado
    if (world_rank == 0) {
        cout << "Proteína traduzida (aminoácidos representados por números):" << endl;
        for (int aminoacid : global_protein) {
            cout << aminoacid << " ";
        }
        cout << endl;
    }

    MPI_Finalize();
    return 0;
}
