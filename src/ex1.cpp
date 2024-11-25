#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>

using namespace std;

void count_nucleotides(const string& filename, int& countA, int& countT, int& countC, int& countG) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Erro ao abrir o arquivo: " << filename << endl;
        return;
    }

    // Contagem local por thread
    vector<int> local_countA, local_countT, local_countC, local_countG;
    #pragma omp parallel
    {
        int localA = 0, localT = 0, localC = 0, localG = 0;

        string line;
        while (getline(file, line)) {
            #pragma omp for schedule(static)
            for (size_t i = 0; i < line.size(); ++i) {
                switch (line[i]) {
                    case 'A': localA++; break;
                    case 'T': localT++; break;
                    case 'C': localC++; break;
                    case 'G': localG++; break;
                }
            }
        }

        #pragma omp critical
        {
            local_countA.push_back(localA);
            local_countT.push_back(localT);
            local_countC.push_back(localC);
            local_countG.push_back(localG);
        }
    }

    // Redução para somar os resultados das threads
    countA = countT = countC = countG = 0;
    for (int count : local_countA) countA += count;
    for (int count : local_countT) countT += count;
    for (int count : local_countC) countC += count;
    for (int count : local_countG) countG += count;

    file.close();
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

    // Lista de arquivos
    vector<string> files;
    for (int i = 1; i <= num_files; ++i) {
        files.push_back("dados/processados/processed_chr" + to_string(i) + ".subst.fa");
    }

    // Dividir arquivos entre processos
    int files_per_process = files.size() / world_size;
    int start_idx = world_rank * files_per_process;
    int end_idx = (world_rank == world_size - 1) ? files.size() : start_idx + files_per_process;

    // Contagem local por processo
    int totalA = 0, totalT = 0, totalC = 0, totalG = 0;
    for (int i = start_idx; i < end_idx; ++i) {
        int countA = 0, countT = 0, countC = 0, countG = 0;
        count_nucleotides(files[i], countA, countT, countC, countG);

        totalA += countA;
        totalT += countT;
        totalC += countC;
        totalG += countG;
    }

    // Redução global para agregar os resultados de todos os processos
    int globalA = 0, globalT = 0, globalC = 0, globalG = 0;
    MPI_Reduce(&totalA, &globalA, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&totalT, &globalT, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&totalC, &globalC, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&totalG, &globalG, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Apenas o processo de rank 0 exibe os resultados
    if (world_rank == 0) {
        cout << "Contagem total de nucleotídeos:" << endl;
        cout << "A: " << globalA << endl;
        cout << "T: " << globalT << endl;
        cout << "C: " << globalC << endl;
        cout << "G: " << globalG << endl;
    }

    MPI_Finalize();
    return 0;
}
