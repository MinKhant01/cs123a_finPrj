import math                                         # for math operations
import argparse                                     # for command line argument parsing   
import numpy as np                                  # creating and manipulating matrices
import time                                         # time tracking for benchmarking
import resource                                     # resource usage for benchmarking
import sys                                          # for safe exit from the program
import os                                           # for file and directory operations
import csv                                          # creating and manipulating CSV files
from Bio.Data.IUPACData import protein_letters      # defining allowed protein letters
from Bio.Align import substitution_matrices         # for loading BLOSUM matrices

# BLOSUM matrices
BLOSUM45 = substitution_matrices.load("BLOSUM45")
BLOSUM50 = substitution_matrices.load("BLOSUM50")
BLOSUM62 = substitution_matrices.load("BLOSUM62")
BLOSUM80 = substitution_matrices.load("BLOSUM80")

def parse_fasta(file_path):
    '''
        input {file_path: str -> path to the FASTA file}
        output {sequence: str -> the protein sequence from the FASTA file}

        This function opens the FASTA file, reads its content, and returns the sequence content,
        and raises an exception if the file is not in a valid FASTA format.
    '''
    try:
        with open(file_path, "r") as file:
            lines = file.read().splitlines()
    except Exception as e:
        raise Exception(f"Error reading file {file_path}: {e}") # raise exception if file not found or unreadable
    
    if not lines or not lines[0].startswith(">"):
        raise Exception(f"File {file_path} does not contain a valid FASTA format.") # raise exception if file is not in FASTA format
    
    sequence = "".join(line.strip() for line in lines[1:] if line.strip() and not line.startswith(">")) # concatenate all lines except the header, blank lines and whitespace
    return sequence

def validate_sequence(sequence):
    '''
        input {sequence: str -> the protein sequence to validate}
        output {bool -> True if the sequence is valid, False otherwise}

        This function validates the protein sequence by checking 
        if it contains only valid amino acid letters as allowed by the IUPAC standard.
    '''
    seq = sequence.upper()                  # convert to uppercase for consistency
    if not seq:                             # check if the sequence is empty
        return False
    allowed_set = set(protein_letters)
    for residue in seq:                     # iterate through each residue in the sequence for validation
        if residue not in allowed_set:
            return False
    return True

def get_substitution_score(res1, res2, matrix):
    '''
        input {res1: str -> first amino acid;
                res2: str -> second amino acid;
                matrix: dict -> substitution matrix (e.g., BLOSUM)}
        output {score: int -> substitution score} 

        This function returns the substitution score for a pair of amino acids
        using the provided substitution matrix. It returns -1 if the pair is not found in the matrix.
    '''
    res1, res2 = res1.upper(), res2.upper()     # convert to uppercase for consistency
    if (res1, res2) in matrix:                  # check if the pair exists in the matrix
        return matrix[(res1, res2)]
    elif (res2, res1) in matrix:                # check if the reverse pair exists in the matrix
        return matrix[(res2, res1)]
    else:
        return -1
    
def get_context_matrix(seq1, seq2, matrix):
    '''
        input {seq1: str -> first protein sequence;
                seq2: str -> second protein sequence;
                matrix: dict -> substitution matrix (e.g., BLOSUM)}
        output {matrix: np.array -> context matrix}

        This function creates a context matrix for the two sequences using the provided substitution matrix.
        The context matrix is a 2D array where each cell (i, j) contains the substitution score
        for the i-th residue of seq1 and the j-th residue of seq2.
    '''
    rows, cols = len(seq1), len(seq2)
    score_matrix = [[0 for j in range(cols)] for i in range(rows)]
    for i in range(rows):
        for j in range(cols):
            score_matrix[i][j] = get_substitution_score(seq1[i], seq2[j], matrix)
    return np.array(score_matrix)

def smith_waterman(seq1, seq2, sub_matrix, gap_penalty=-5):
    '''
        input {seq1: str -> first protein sequence;
                seq2: str -> second protein sequence;
                sub_matrix: dict -> substitution matrix (e.g., BLOSUM);
                gap_penalty: int -> gap penalty (default: -5)}
        output {align1: str -> aligned sequence 1;
                align2: str -> aligned sequence 2;
                max_score: int -> maximum alignment score;
                H: np.array -> DP matrix;
                max_pos: tuple -> position of the maximum score in the DP matrix;
                start_pos: tuple -> starting position of the alignment in the DP matrix}

        This function performs the Smith-Waterman local alignment algorithm
        using the provided substitution matrix and gap penalty. It returns 
            the aligned sequences,
            the maximum alignment score, 
            the DP matrix, 
            the position of the maximum score in the DP matrix, and
            the starting position of the alignment in the DP matrix.
    '''
    rows, cols = len(seq1), len(seq2)
    max_score = 0
    max_pos = (0, 0)

    # initialize DP matrix
    H = np.zeros((rows+1, cols+1), dtype=int)
    
    # fill in DP matrix
    for i in range(1, rows+1):
        for j in range(1, cols+1):
            sub_score = get_substitution_score(seq1[i-1], seq2[j-1], sub_matrix)
            diag = H[i-1][j-1] + sub_score
            up = H[i-1][j] + gap_penalty
            left = H[i][j-1] + gap_penalty

            H[i][j] = max(0, diag, up, left)
            if H[i][j] > max_score:
                max_score = H[i][j]
                max_pos = (i, j)

    # traceback to find the optimal local alignment
    align1, align2 = "", ""
    i, j = max_pos
    while i > 0 and j > 0 and H[i][j] > 0:
        score_current = H[i][j]
        score_diag = H[i-1][j-1]
        score_up = H[i-1][j]
        score_left = H[i][j-1]
        sub_score = get_substitution_score(seq1[i-1], seq2[j-1], sub_matrix)
        if score_current == score_diag + sub_score:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif score_current == score_up + gap_penalty:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        elif score_current == score_left + gap_penalty:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
        else:
            break
    return align1, align2, max_score, H, max_pos, (i, j)

def main():
    '''
        This function parses command line arguments, checks if FASTA files exist, parses them to extract sequences,,
        validates the sequences, performs Smith-Waterman alignment, creates a folder for each command, 
        writes output files to the folder in the project directory. The output files include:
            - context matrix CSV file
            - DP matrix CSV file
            - metrics text file
        The metrics text file contains:
            - alignment score
            - aligned sequences
            - percentage identity
            - coverage
            - gap metrics
            - estimated E-value
            - benchmark results (elapsed time, CPU time, memory usage)
    '''

    # parse command line arguments; takes in two FASTA files, a BLOSUM matrix, and a gap penalty
    parser = argparse.ArgumentParser(description="Smith-Waterman local alignment using a BLOSUM matrix.")
    parser.add_argument("file1", type=str, help="FASTA file containing protein sequence 1")
    parser.add_argument("file2", type=str, help="FASTA file containing protein sequence 2")
    parser.add_argument("--matrix", type=str, default="BLOSUM62", choices=["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "ALL"], help="BLOSUM matrix family to use (default: BLOSUM62)")
    parser.add_argument("--gap_penalty", type=int, default=-5, help="Gap penalty (default: -5)")
    parser.add_argument("--ungapped", action="store_true", default=False , help="Use ungapped alignment (default: False)")
    args = parser.parse_args()

    # check if the input files exist
    if not os.path.isfile(args.file1):
        print(f"Error: File {args.file1} does not exist.")
        sys.exit(1)
    if not os.path.isfile(args.file2):
        print(f"Error: File {args.file2} does not exist.")
        sys.exit(1)
    
    # validate the FASTA files
    if not args.file1.endswith(".fasta") and not args.file1.endswith(".fa"):
        print(f"Error: File {args.file1} is not a valid FASTA file.")
        sys.exit(1)
    if not args.file2.endswith(".fasta") and not args.file2.endswith(".fa"):
        print(f"Error: File {args.file2} is not a valid FASTA file.")
        sys.exit(1)

    # parse the FASTA files and extract the sequences
    try:                                            
        seq1 = parse_fasta(args.file1)
    except Exception as e:
        print(f"Error parsing FASTA file {args.file1}: {e}")
        sys.exit(1)
    try:
        seq2 = parse_fasta(args.file2)
    except Exception as e:
        print(f"Error parsing FASTA file {args.file2}: {e}")
        sys.exit(1)
    
    # validate the sequences
    if not validate_sequence(seq1):
        print(f"Error: Sequence 1 - {seq1} - is not a valid protein sequence.")
        return
    if not validate_sequence(seq2):
        print(f"Error: Sequence 2 - {seq2} - is not a valid protein sequence.")
        return
    
    gap_penalty = args.gap_penalty
    matrix_input = args.matrix.upper()

    # Check if the user wants to run all matrices or a specific one
    if matrix_input not in ["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "ALL"]:
        print(f"Error: {matrix_input} is currently not a supported BLOSUM matrix.")
        sys.exit(1)

    # run all matrices if ALL is specified by the user, or run the specified matrix otherwise
    if matrix_input == "ALL":
        matrices_to_run = ["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80"]
    else:
        matrices_to_run = [matrix_input]

    # interate through the matrices/matrix to perform Smith-Waterman alignment
    for matrix_input in matrices_to_run:
        match matrix_input:
            case "BLOSUM50":
                sub_matrix = BLOSUM50
            case "BLOSUM62":
                sub_matrix = BLOSUM62
            case "BLOSUM80":
                sub_matrix = BLOSUM80
            case "BLOSUM45":
                sub_matrix = BLOSUM45
            case _:
                print(f"Error: {matrix_input} is currently not a supported BLOSUM matrix.")
                continue
    
        context_matrix = get_context_matrix(seq1, seq2, sub_matrix)

        # extract base names from file paths
        file1_base = os.path.splitext(os.path.basename(args.file1))[0]
        file2_base = os.path.splitext(os.path.basename(args.file2))[0]

        # sort bases alphabetically for consistent ordering
        bases = sorted([file1_base, file2_base])
        output_base = f"{bases[0]}-{bases[1]}"
        output_dir = f"{output_base}_{matrix_input}_output"
        os.makedirs(output_dir, exist_ok=True)

        # write out the context matrix as a CSV file
        context_csv_filename = os.path.join(output_dir, f"{output_base}_{matrix_input}_contextMatrix.csv")
        with open(context_csv_filename, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            for row in context_matrix:
                writer.writerow(row)

        # record benchmark metrics
        start_time = time.time()
        align1, align2, score, dp_matrix, max_pos, start_pos = smith_waterman(seq1, seq2, sub_matrix, gap_penalty)
        end_time = time.time()
        elapsed_time = end_time - start_time
        usage = resource.getrusage(resource.RUSAGE_SELF)

        # percentage Identity: count of positions with identical amino acids (ignoring gaps) relative to alignment length.
        identity_matches = sum(1 for a, b in zip(align1, align2) if a == b and a != "-")
        percent_identity = (identity_matches / len(align1)) * 100 if align1 else 0

        # coverage: proportion of each sequence that is aligned.
        # since the traceback stops at the cell where H is 0, the aligned segment spans from start_pos to max_pos.
        coverage_seq1 = ((max_pos[0] - start_pos[0]) / len(seq1)) * 100
        coverage_seq2 = ((max_pos[1] - start_pos[1]) / len(seq2)) * 100

        # gap metrics: count and percentage of gaps in each aligned sequence.
        gaps_seq1 = align1.count("-")
        gaps_seq2 = align2.count("-")
        gap_percentage_seq1 = (gaps_seq1 / len(align1)) * 100 if align1 else 0
        gap_percentage_seq2 = (gaps_seq2 / len(align2)) * 100 if align2 else 0

        # e-value estimation
        # e-value = K * m * n * exp(-lambda_val * score)
        if args.ungapped:
            K = 0.041           # default ungapped value
            lambda_val = 0.267  # default ungapped value
        else:
            K = 0.128           # default gapped value
            lambda_val = 0.311  # default gapped value
        m = len(seq1)
        n = len(seq2)
        e_value = K * m * n * math.exp(-lambda_val * score)

        # write out the DP matrix as a CSV file
        dp_csv_filename = os.path.join(output_dir, f"{output_base}_{matrix_input}_DPMatrix.csv")
        with open(dp_csv_filename, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            for row in dp_matrix:
                writer.writerow(row)

        # write out the metrics (alignment results and benchmark) as a txt file
        metrics_filename = os.path.join(output_dir, f"{output_base}_{matrix_input}_metrics.txt")
        with open(metrics_filename, "w") as f:
            # alignment results
            f.write(f"Alignment Score: {score}\n")
            f.write(f"Aligned Sequence 1: {align1}\n")
            f.write(f"Aligned Sequence 2: {align2}\n")
            f.write(f"Percentage Identity: {percent_identity:.2f}%\n")
            f.write(f"Coverage - Sequence 1: {coverage_seq1:.2f}%, Sequence 2: {coverage_seq2:.2f}%\n")
            f.write(f"Gap Metrics - Sequence 1: {gaps_seq1} gaps ({gap_percentage_seq1:.2f}%), ")
            f.write(f"Sequence 2: {gaps_seq2} gaps ({gap_percentage_seq2:.2f}%)\n")
            f.write(f"Estimated E-value: {e_value:.4e}\n\n")
            # benchmark results
            f.write("Benchmark Results:\n")
            f.write(f"Elapsed (wall-clock) Time: {elapsed_time:.4f} seconds\n")
            f.write(f"CPU Time: User {usage.ru_utime:.4f} sec, System {usage.ru_stime:.4f} sec\n")
            f.write(f"Maximum Memory Usage: {usage.ru_maxrss} KB\n")
        
        # prints an ack message to the console
        print(f"Files generated:\n {context_csv_filename}\n {dp_csv_filename}\n {metrics_filename}")

if __name__ == "__main__":
    main()


    