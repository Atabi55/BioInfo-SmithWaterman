# Question 1: Needleman-Wunsch Global Alignment with match = 2, mismatch = -1, and g = -2

import numpy as np

# Define the sequences to be aligned
seq1 = "AATCGTTCGAC"
seq2 = "AACGTTTTCGGCA"

# Define scoring scheme
match_score = 2
mismatch_score = -1
g = -2  # gap penalty

# Initialize the scoring matrix with zeros
rows = len(seq2) + 1
cols = len(seq1) + 1
matrix = np.zeros((rows, cols), dtype=int)

# Fill the first row and column with cumulative gap penalties
for i in range(1, rows):
    matrix[i][0] = matrix[i-1][0] + g
for j in range(1, cols):
    matrix[0][j] = matrix[0][j-1] + g

# Fill the matrix using the Needleman-Wunsch recurrence relation
for i in range(1, rows):
    for j in range(1, cols):
        if seq2[i - 1] == seq1[j - 1]:
            score = match_score  # match
        else:
            score = mismatch_score  # mismatch

        diagonal = matrix[i - 1][j - 1] + score  # align characters
        up = matrix[i - 1][j] + g                # gap in seq1
        left = matrix[i][j - 1] + g              # gap in seq2

        matrix[i][j] = max(diagonal, up, left)  # choose the best move

# Traceback to find the optimal alignment
aligned_seq1 = ""
aligned_seq2 = ""
i = len(seq2)
j = len(seq1)

# Follow the path backwards from bottom-right to top-left
while i > 0 or j > 0:
    current = matrix[i][j]

    if i > 0 and j > 0:
        if seq2[i - 1] == seq1[j - 1]:
            score = match_score
        else:
            score = mismatch_score

        if current == matrix[i - 1][j - 1] + score:
            aligned_seq1 = seq1[j - 1] + aligned_seq1
            aligned_seq2 = seq2[i - 1] + aligned_seq2
            i -= 1
            j -= 1
            continue

    if i > 0 and current == matrix[i - 1][j] + g:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[i - 1] + aligned_seq2
        i -= 1
        continue

    if j > 0:
        aligned_seq1 = seq1[j - 1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        j -= 1

# Output the result of the alignment
print("Aligned Sequence 1:", aligned_seq1)
print("Aligned Sequence 2:", aligned_seq2)
print("Alignment Score:", matrix[-1][-1])