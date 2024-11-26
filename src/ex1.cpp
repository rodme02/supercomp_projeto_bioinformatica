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

    // Read entire file content into a single string
    string file_content((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
    file.close();

    // Initialize counts
    countA = countT = countC = countG = 0;

    // Parallelize counting
    #pragma omp parallel for schedule(static) reduction(+:countA, countT, countC, countG)
    for (size_t i = 0; i < file_content.size(); ++i) {
        switch (file_content[i]) {
            case 'A': countA++; break;
            case 'T': countT++; break;
            case 'C': countC++; break;
            case 'G': countG++; break;
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Default number of files to process
    int num_files = 22;
    if (argc > 1) {
        num_files = stoi(argv[1]);
    }

    // List of files
    vector<string> files;
    for (int i = 1; i <= num_files; ++i) {
        files.push_back("dados/processados/processed_chr" + to_string(i) + ".subst.fa");
    }

    // Divide files among processes
    int files_per_process = files.size() / world_size;
    int start_idx = world_rank * files_per_process;
    int end_idx = (world_rank == world_size - 1) ? files.size() : start_idx + files_per_process;

    // Local nucleotide counts for each process
    int totalA = 0, totalT = 0, totalC = 0, totalG = 0;
    for (int i = start_idx; i < end_idx; ++i) {
        int countA = 0, countT = 0, countC = 0, countG = 0;
        count_nucleotides(files[i], countA, countT, countC, countG);

        // Print the counts for each file
        cout << "Processo " << world_rank << " processou arquivo: " << files[i] << endl;
        cout << "Contagens no arquivo: A: " << countA
             << ", T: " << countT
             << ", C: " << countC
             << ", G: " << countG << endl;

        totalA += countA;
        totalT += countT;
        totalC += countC;
        totalG += countG;
    }

    // Global reduction to aggregate results from all processes
    int globalA = 0, globalT = 0, globalC = 0, globalG = 0;
    MPI_Reduce(&totalA, &globalA, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&totalT, &globalT, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&totalC, &globalC, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&totalG, &globalG, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Only rank 0 displays global results
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
