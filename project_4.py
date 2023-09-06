"""
This module contains functions for measuring the similarity between two sequences of characters.
It uses dynamic programming for computing an alignment table based on the values of their scoring matrix.
"""


def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    Takes as input a set of characters (alphabet) and three score values. 
    Returns a dictionary of dictionaries whose entries are indexed by pairs of characters in (alphabet) plus '-'.
    Every cell in the matrix is scored based on its position.
    The score for any entry indexed by one or more dashes is always (dash_score). 
    """
    entries = alphabet.copy()
    entries.add("-")
    scoring_matrix = {char_1 : {char_2 : 0} for char_1 in entries for char_2 in entries}
    for matrix_row in entries:
        for matrix_col in entries:
            if matrix_row == "-" or matrix_col == "-":
                score = dash_score
            elif matrix_row == matrix_col:
                score = diag_score
            else:
                score = off_diag_score

            scoring_matrix[matrix_row][matrix_col] = score

    return scoring_matrix


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    Takes as input two sequences (seq_x) and (seq_y) whose elements share a common alphabet with the (scoring matrix). 
    Returns a list of lists representing the dynamic programming table for these sequences. 
    If (global_flag) is "True", each entry of the table is computed for global allignment.
    Otherwise the local allignment computation method is used.
    """
    rows = range(len(seq_x) + 1)
    cols = range(len(seq_y) + 1)

    align_matrix = []
    for row in rows:
        align_matrix.append([])
        for col in cols:
            if row == 0 and col == 0:
                score = 0
            elif row == 0 and col > 0:
                score = align_matrix[0][col - 1] + scoring_matrix["-"][seq_y[col - 1]]
            elif row > 0 and col == 0:
                score = align_matrix[row - 1][0] + scoring_matrix[seq_x[row - 1]]["-"]
            else:
                up = align_matrix[row - 1][col] + scoring_matrix[seq_x[row - 1]]["-"]
                left = align_matrix[row][col - 1] + scoring_matrix["-"][seq_y[col - 1]]
                diagonal = align_matrix[row - 1][col - 1] + scoring_matrix[seq_x[row - 1]][seq_y[col - 1]]
                score = max(up, left, diagonal)
            if not global_flag:
                if score < 0:
                    score = 0
            align_matrix[row].append(score)

    return align_matrix


def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Takes as input two sequences (seq_x) and (seq_y) whose elements share 
    a common alphabet with the (scoring matrix). 
    Computes a global alignment of the sequences using the global (alignment matrix).
    """
    score = alignment_matrix[-1][-1]
    row, col = len(seq_x), len(seq_y)
    alig_x, alig_y = "", ""

    while row != 0 and col != 0:
        if alignment_matrix[row][col] == alignment_matrix[row - 1][col - 1] + scoring_matrix[seq_x[row - 1]][seq_y[col - 1]]:
            alig_x = seq_x[row - 1] + alig_x
            alig_y = seq_y[col - 1] + alig_y
            row -= 1
            col -= 1
        elif alignment_matrix[row][col] == alignment_matrix[row][col - 1] + scoring_matrix["-"][seq_y[col - 1]]:
            alig_x = "-" + alig_x
            alig_y = seq_y[col - 1] + alig_y
            col -= 1
        else:
            alig_x = seq_x[row - 1] + alig_x
            alig_y = "-" + alig_y
            row -= 1

    while row != 0:
            alig_x = seq_x[row - 1] + alig_x
            alig_y = "-" + alig_y
            row -= 1
    while col != 0:
            alig_x = "-" + alig_x
            alig_y = seq_y[col - 1] + alig_y
            col -= 1

    return score, alig_x, alig_y


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Takes as input two sequences (seq_x) and (seq_y) whose elements share 
    a common alphabet with the (scoring matrix). 
    Computes a local optimal alignment of the sequences using the local (alignment matrix).
    """
    alig_x, alig_y = "", ""
    max_score = float("-inf")

    for idx_1 in range(len(seq_x) + 1):
        for idx_2 in range(len(seq_y) + 1):
            if alignment_matrix[idx_1][idx_2] > max_score:
                max_score = alignment_matrix[idx_1][idx_2]
                row = idx_1
                col = idx_2

    while alignment_matrix[row][col] != 0:
        if alignment_matrix[row][col] == alignment_matrix[row - 1][col - 1] + scoring_matrix[seq_x[row - 1]][seq_y[col - 1]]:
            alig_x = seq_x[row - 1] + alig_x
            alig_y = seq_y[col - 1] + alig_y
            row -= 1
            col -= 1
        elif alignment_matrix[row][col] == alignment_matrix[row][col - 1] + scoring_matrix["-"][seq_y[col - 1]]:
            alig_x = "-" + alig_x
            alig_y = seq_y[col - 1] + alig_y
            col -= 1
        else:
            alig_x = seq_x[row - 1] + alig_x
            alig_y = "-" + alig_y
            row -= 1

    return max_score, alig_x, alig_y
