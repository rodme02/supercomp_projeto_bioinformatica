#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <mpi.h>
#include <omp.h>

using namespace std;

// Map codons to amino acids (represented by numbers)
unordered_map<string, int> codon_to_aminoacid = {
    {"CCA", 1}, {"CCG", 1}, {"CCU", 1}, {"CCC", 1}, // Proline
    {"UCU", 2}, {"UCA", 2}, {"UCG", 2}, {"UCC", 2}, // Serine
    {"CAG", 3}, {"CAA", 3},                         // Glutamine
    {"ACA", 4}, {"ACC", 4}, {"ACU", 4}, {"ACG", 4}, // Threonine
    {"AUG", 5},                                     // Methionine (Start)
    {"UGC", 6}, {"UGU", 6},                         // Cysteine
    {"GUG", 7}, {"GUA", 7}, {"GUC", 7}, {"GUU", 7}, // Valine
    {"UGA", -1}                                     // Stop codon
};

// Translate RNA sequence into amino acids
vector<int> translate_rna(const string& rna_sequence) {
    vector<int> protein;
    for (size_t i = 0; i + 3 <= rna_sequence.size(); i += 3) { // Ensure at least 3 characters remain
        string codon = rna_sequence.substr(i, 3);
        if (codon_to_aminoacid.find(codon) != codon_to_aminoacid.end()) {
            int aminoacid = codon_to_aminoacid[codon];
            if (aminoacid == -1) break; // Stop codon
            protein.push_back(aminoacid);
        }
    }
    return protein;
}


// Process an RNA file
vector<int> process_file(const string& file_path) {
    ifstream file(file_path);
    if (!file.is_open()) {
        cerr << "Error opening file: " << file_path << endl;
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

    int num_files = 22;
    if (argc > 1) {
        num_files = stoi(argv[1]);
    }

    vector<string> rna_files;
    for (int i = 1; i <= num_files; ++i) {
        rna_files.push_back("dados/rna/rna_chr" + to_string(i) + ".subst.fa");
    }

    int files_per_process = rna_files.size() / world_size;
    int start_idx = world_rank * files_per_process;
    int end_idx = (world_rank == world_size - 1) ? rna_files.size() : start_idx + files_per_process;

    vector<int> local_protein;
    for (int i = start_idx; i < end_idx; ++i) {
        cout << "Process " << world_rank << " processing file: " << rna_files[i] << endl;
        vector<int> protein = process_file(rna_files[i]);
        local_protein.insert(local_protein.end(), protein.begin(), protein.end());
    }

    // Step 1: Gather the sizes of local_protein from all processes
    int local_size = local_protein.size();
    vector<int> all_sizes;
    if (world_rank == 0) {
        all_sizes.resize(world_size);
    }
    MPI_Gather(&local_size, 1, MPI_INT,
            world_rank == 0 ? all_sizes.data() : NULL, 1, MPI_INT,
            0, MPI_COMM_WORLD);

    // Step 2: Calculate displacement for MPI_Gatherv
    vector<int> displacements;
    int total_size = 0;
    if (world_rank == 0) {
        displacements.resize(world_size, 0);
        for (int i = 0; i < world_size; ++i) {
            displacements[i] = total_size;
            total_size += all_sizes[i];
        }
    }

    // Step 3: Gather all proteins into global_protein
    vector<int> global_protein;
    if (world_rank == 0) {
        global_protein.resize(total_size);
    }
    MPI_Gatherv(local_protein.data(), local_size, MPI_INT,
                world_rank == 0 ? global_protein.data() : NULL,
                world_rank == 0 ? all_sizes.data() : NULL,
                world_rank == 0 ? displacements.data() : NULL,
                MPI_INT, 0, MPI_COMM_WORLD);

    // Step 4: Print results on the master process
    if (world_rank == 0) {
        cout << "Translated protein (amino acids represented by numbers):" << endl;
        cout << "Total amino acids: " << global_protein.size() << endl;

        // Limit the number of amino acids to print
        size_t max_print = 100;
        for (size_t i = 0; i < min(global_protein.size(), max_print); ++i) {
            cout << global_protein[i] << " ";
        }
        if (global_protein.size() > max_print) {
            cout << "... (output truncated)" << endl;
        } else {
            cout << endl;
        }
        cout << "Finished printing amino acids." << endl;

        // Verify data integrity
        bool invalid_data = false;
        for (int aminoacid : global_protein) {
            if (aminoacid < -1 || aminoacid > 7) {
                invalid_data = true;
                cout << "Invalid amino acid detected: " << aminoacid << endl;
                break;
            }
        }
        if (!invalid_data) {
            cout << "All amino acids are valid." << endl;
        }
    }

    cout << "Process " << world_rank << " reached the end of the program." << endl;
    cout.flush();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
