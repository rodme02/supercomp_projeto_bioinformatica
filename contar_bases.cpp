#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cctype>

// Função para contar as bases em uma sequência de DNA (excluindo caracteres 'N')
std::map<char, int> contar_bases(const std::string& dna_seq) {
    std::map<char, int> contagem = {{'A', 0}, {'T', 0}, {'C', 0}, {'G', 0}};
    for (char base : dna_seq) {
        base = std::toupper(base); // Converte para maiúsculas
        if (contagem.find(base) != contagem.end()) {
            contagem[base]++;
        }
    }
    return contagem;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Nome do arquivo de DNA
    std::string filename = "dados/chr1.subst.fa"; // Modifique para cada cromossomo

    std::string dna_seq;
    if (rank == 0) {
        // Processo mestre lê o arquivo de DNA e filtra as bases
        std::ifstream file(filename);
        if (file) {
            std::string line;
            while (getline(file, line)) {
                // Ignora a linha de descrição e linhas contendo apenas 'N'
                if (line[0] == '>' || line.find_first_not_of("N") == std::string::npos) {
                    continue;
                }
                // Converte para maiúsculas e adiciona à sequência
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                dna_seq += line;
            }
        }
    }

    // Enviar a sequência filtrada para todos os processos
    int seq_size = dna_seq.size() / size;
    std::string local_seq(seq_size, ' ');
    MPI_Scatter(dna_seq.c_str(), seq_size, MPI_CHAR, &local_seq[0], seq_size, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Cada processo conta suas bases
    std::map<char, int> local_count = contar_bases(local_seq);

    // Agregar os resultados no processo mestre
    std::map<char, int> global_count;
    for (auto& pair : local_count) {
        int total_count;
        MPI_Reduce(&pair.second, &total_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            global_count[pair.first] = total_count;
        }
    }

    // Exibir resultados no processo mestre
    if (rank == 0) {
        std::cout << "Contagem de bases:\n";
        for (const auto& pair : global_count) {
            std::cout << pair.first << ": " << pair.second << '\n';
        }
    }

    MPI_Finalize();
    return 0;
}
