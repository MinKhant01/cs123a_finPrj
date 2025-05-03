import argparse
from Bio.Data.IUPACData import protein_letters
from Bio.Align import substitution_matrices
import numpy as np

# BLOSUM matrices
BLOSUM45 = substitution_matrices.load("BLOSUM45")
BLOSUM50 = substitution_matrices.load("BLOSUM50")
BLOSUM62 = substitution_matrices.load("BLOSUM62")
BLOSUM80 = substitution_matrices.load("BLOSUM80")

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
    
def get_substitution_matrix(seq1, seq2, matrix):
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
    return align1, align2, max_score, H

def main():
    parser = argparse.ArgumentParser(description="Smith-Waterman local alignment using a BLOSUM matrix.")
    parser.add_argument("seq1", type=str, help="protein sequence 1")
    parser.add_argument("seq2", type=str, help="protein sequence 2")
    parser.add_argument("--matrix", type=str, default="BLOSUM62", choices=["BLOSUM30", "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80"], help="BLOSUM matrix family to use (default: BLOSUM62)")
    parser.add_argument("--gap_penalty", type=int, default=-5, help="Gap penalty (default: -5)")
    args = parser.parse_args()

    if not validate_sequence(args.seq1):
        print(f"Error: Sequence 1 - {args.seq1} - is not a valid protein sequence.")
        return
    if not validate_sequence(args.seq2):
        print(f"Error: Sequence 2 - {args.seq2} - is not a valid protein sequence.")
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
    
    context_matrix = get_substitution_matrix(args.seq1, args.seq2, sub_matrix)
    print("Context Scoring Matrix:")
    for row in context_matrix:
        print(row)


    align1, align2, score, dp_matrix = smith_waterman(args.seq1, args.seq2, sub_matrix, gap_penalty)

    print(f"Alignment Score: {score}")
    print(f"Aligned Sequence 1: {align1}")
    print(f"Aligned Sequence 2: {align2}")
    print("DP Matrix:")
    for row in dp_matrix:
        print(row)

if __name__ == "__main__":
    main()


    