import math
import argparse
import numpy as np
import time
import resource
import sys
import os
import csv
from Bio.Data.IUPACData import protein_letters
from Bio.Align import substitution_matrices

# BLOSUM matrices
BLOSUM45 = substitution_matrices.load("BLOSUM45")
BLOSUM50 = substitution_matrices.load("BLOSUM50")
BLOSUM62 = substitution_matrices.load("BLOSUM62")
BLOSUM80 = substitution_matrices.load("BLOSUM80")

def parse_fasta(file_path):
    try:
        with open(file_path, "r") as file:
            lines = file.read().splitlines()
    except Exception as e:
        raise Exception(f"Error reading file {file_path}: {e}")
    
    if not lines or not lines[0].startswith(">"):
        raise Exception(f"File {file_path} does not contain a valid FASTA format.")
    
    sequence = "".join(line.strip() for line in lines[1:] if line.strip() and not line.startswith(">"))
    return sequence

def validate_sequence(sequence):
    seq = sequence.upper()
    allowed_set = set(protein_letters)
    for residue in seq:
        if residue not in allowed_set:
            return False
    return True

def get_substitution_score(res1, res2, matrix):
    res1, res2 = res1.upper(), res2.upper()
    if (res1, res2) in matrix:
        return matrix[(res1, res2)]
    elif (res2, res1) in matrix:
        return matrix[(res2, res1)]
    else:
        return -1
    
def get_context_matrix(seq1, seq2, matrix):
    rows, cols = len(seq1), len(seq2)
    score_matrix = [[0 for j in range(cols)] for i in range(rows)]
    for i in range(rows):
        for j in range(cols):
            score_matrix[i][j] = get_substitution_score(seq1[i], seq2[j], matrix)
    return np.array(score_matrix)

def smith_waterman(seq1, seq2, sub_matrix, gap_penalty=-5):
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
    parser = argparse.ArgumentParser(description="Smith-Waterman local alignment using a BLOSUM matrix.")
    parser.add_argument("file1", type=str, help="FASTA file containing protein sequence 1")
    parser.add_argument("file2", type=str, help="FASTA file containing protein sequence 2")
    parser.add_argument("log_file", type=str, help="log file to record the output")
    parser.add_argument("--matrix", type=str, default="BLOSUM62", choices=["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80"], help="BLOSUM matrix family to use (default: BLOSUM62)")
    parser.add_argument("--gap_penalty", type=int, default=-5, help="Gap penalty (default: -5)")
    parser.add_argument("--ungapped", action="store_true", default=False , help="Use ungapped alignment (default: False)")
    args = parser.parse_args()

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
    

    if not validate_sequence(seq1):
        print(f"Error: Sequence 1 - {seq1} - is not a valid protein sequence.")
        return
    if not validate_sequence(seq2):
        print(f"Error: Sequence 2 - {seq2} - is not a valid protein sequence.")
        return

    matrix_input = args.matrix.upper()
    gap_penalty = args.gap_penalty

    match matrix_input:
        case "BLOSUM45":
            sub_matrix = BLOSUM45
        case "BLOSUM50":
            sub_matrix = BLOSUM50
        case "BLOSUM62":
            sub_matrix = BLOSUM62
        case "BLOSUM80":
            sub_matrix = BLOSUM80
        case _:
            print(f"Error: {matrix_input} is currently not a supported BLOSUM matrix.")
            sys.exit(1)
    
    context_matrix = get_context_matrix(seq1, seq2, sub_matrix)

    # Derive an output base name from file1 (without extension)
    output_base = os.path.splitext(os.path.basename(args.log_file))[0]
    output_dir = output_base + "_output"
    os.makedirs(output_dir, exist_ok=True)

    context_csv_filename = os.path.join(output_dir, f"{output_base}_contextMatrix.csv")
    with open(context_csv_filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        for row in context_matrix:
            writer.writerow(row)

    # Record the start wall-clock time.
    start_time = time.time()
    align1, align2, score, dp_matrix, max_pos, start_pos = smith_waterman(seq1, seq2, sub_matrix, gap_penalty)
    end_time = time.time()
    elapsed_time = end_time - start_time
    usage = resource.getrusage(resource.RUSAGE_SELF)

    # Percentage Identity: count of positions with identical amino acids (ignoring gaps) relative to alignment length.
    identity_matches = sum(1 for a, b in zip(align1, align2) if a == b and a != "-")
    percent_identity = (identity_matches / len(align1)) * 100 if align1 else 0

    # Coverage: proportion of each sequence that is aligned.
    # Since the traceback stops at the cell where H == 0, the aligned segment spans from start_pos to max_pos.
    coverage_seq1 = ((max_pos[0] - start_pos[0]) / len(seq1)) * 100
    coverage_seq2 = ((max_pos[1] - start_pos[1]) / len(seq2)) * 100

    # Gap Metrics: count and percentage of gaps in each aligned sequence.
    gaps_seq1 = align1.count("-")
    gaps_seq2 = align2.count("-")
    gap_percentage_seq1 = (gaps_seq1 / len(align1)) * 100 if align1 else 0
    gap_percentage_seq2 = (gaps_seq2 / len(align2)) * 100 if align2 else 0

    # Estimate the E-value
    # E-value = K * m * n * exp(-λ * score)
    if args.ungapped:
        K = 0.041           # default ungapped value
        lambda_val = 0.267  # default ungapped value
    else:
        K = 0.128           # default gapped value
        lambda_val = 0.311  # default gapped value
    m = len(seq1)
    n = len(seq2)
    e_value = K * m * n * math.exp(-lambda_val * score)

   # Write out the DP matrix as a CSV file.
    dp_csv_filename = os.path.join(output_dir, f"{output_base}_DPMatrix.csv")
    with open(dp_csv_filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        for row in dp_matrix:
            writer.writerow(row)

    # Write out the metrics (alignment results and benchmark).
    metrics_filename = os.path.join(output_dir, f"{output_base}_metrics.txt")
    with open(metrics_filename, "w") as f:
        # Alignment results
        f.write(f"Alignment Score: {score}\n")
        f.write(f"Aligned Sequence 1: {align1}\n")
        f.write(f"Aligned Sequence 2: {align2}\n")
        f.write(f"Percentage Identity: {percent_identity:.2f}%\n")
        f.write(f"Coverage - Sequence 1: {coverage_seq1:.2f}%, Sequence 2: {coverage_seq2:.2f}%\n")
        f.write(f"Gap Metrics - Sequence 1: {gaps_seq1} gaps ({gap_percentage_seq1:.2f}%), ")
        f.write(f"Sequence 2: {gaps_seq2} gaps ({gap_percentage_seq2:.2f}%)\n")
        f.write(f"Estimated E-value: {e_value:.4e}\n\n")
        # Benchmark results
        f.write("Benchmark Results:\n")
        f.write(f"Elapsed (wall-clock) Time: {elapsed_time:.4f} seconds\n")
        f.write(f"CPU Time: User {usage.ru_utime:.4f} sec, System {usage.ru_stime:.4f} sec\n")
        f.write(f"Maximum Memory Usage: {usage.ru_maxrss} KB\n")
    
    # Also print a confirmation message to the terminal.
    print(f"Files generated:\n {context_csv_filename}\n {dp_csv_filename}\n {metrics_filename}")

if __name__ == "__main__":
    main()


    