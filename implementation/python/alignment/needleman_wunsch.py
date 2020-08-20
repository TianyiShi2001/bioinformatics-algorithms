import sys
import numpy as np

MATCH = 2
MISMATCH = -1
GAP = -1
GAP_CHAR = "-"

# -----------------------------jjj-----------------------------//

SUBSTITUTION_MATRIX = {
    "A": {"A": MATCH, "C": MISMATCH, "G": MISMATCH, "T": MISMATCH},
    "C": {"A": MISMATCH, "C": MATCH, "G": MISMATCH, "T": MISMATCH},
    "G": {"A": MISMATCH, "C": MISMATCH, "G": MATCH, "T": MISMATCH},
    "T": {"A": MISMATCH, "C": MISMATCH, "G": MISMATCH, "T": MATCH},
}

# Set the constants that represent the three directions, which are used in the traceback matrix
UP = 1
LEFT = 2
DIAG = 3


def init_matrices(nrow, ncol):
    score_matrix = np.zeros((nrow, ncol), np.int32)
    traceback_matrix = np.zeros((nrow, ncol), np.int32)
    for i in range(1, nrow):
        score_matrix[i][0] = score_matrix[i - 1][0] + GAP
        traceback_matrix[i][0] = UP
    for j in range(1, ncol):
        score_matrix[0][j] = score_matrix[0][j - 1] + GAP
        traceback_matrix[0][j] = LEFT
    return score_matrix, traceback_matrix


def compute_score_and_traceback_matrices(s1, s2, substitution_matrix):
    nrow, ncol = len(s1) + 1, len(s2) + 1
    score_matrix, traceback_matrix = init_matrices(nrow, ncol)
    for i in range(1, nrow):
        for j in range(1, ncol):
            up = score_matrix[i - 1][j] + GAP
            left = score_matrix[i][j - 1] + GAP
            diag = (
                score_matrix[i - 1][j - 1] + substitution_matrix[s1[i - 1]][s2[j - 1]]
            )
            max_score, direction = compute_max_score_and_direction(up, left, diag)
            score_matrix[i][j] = max_score
            traceback_matrix[i][j] = direction
    return score_matrix, traceback_matrix


def traceback(traceback_matrix, s1, s2, gap_char):
    s1, s2 = list(s1), list(s2)
    i, j = len(s1), len(s2)
    aln1, aln2 = [], []
    while i > 0 and j > 0:
        direction = traceback_matrix[i][j]
        if direction == UP:
            aln1.append(s1.pop())
            aln2.append(gap_char)
            i -= 1
        elif direction == LEFT:
            aln1.append(gap_char)
            aln2.append(s2.pop())
            j -= 1
        else:
            aln1.append(s1.pop())
            aln2.append(s2.pop())
            i -= 1
            j -= 1
    return ("".join(reversed(aln1)), "".join(reversed(aln2)))


def needleman_wunsch(
    s1, s2, gap_char=GAP_CHAR, substitution_matrix=SUBSTITUTION_MATRIX
):
    score_matrix, traceback_matrix = compute_score_and_traceback_matrices(
        s1, s2, substitution_matrix
    )
    max_score = score_matrix[-1, -1]
    aln1, aln2 = traceback(traceback_matrix, s1, s2, gap_char)
    return (max_score, aln1, aln2)


def compute_max_score_and_direction(up, left, diag):
    max_score = max(up, left, diag)
    if max_score == up:
        direction = UP
    elif max_score == left:
        direction = LEFT
    else:
        direction = DIAG
    return max_score, direction


if __name__ == "__main__":
    s1, s2 = sys.argv[1:3]
    print(needleman_wunsch(s1, s2))

