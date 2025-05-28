




def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2, use_blosum62=False):
    #Inputs: 2 DNA sequences, match, mismatch, gap scores, use_blosum62 boolean
    #Default values: match=2, mismatch=-1, gap=-2, use_blosum62=False
    #Output: List of all local alignments with inserted gaps


    #matrix: 2-dimensional list to hold score of each cell
    #seq1 in vertical axis, seq2 in horizontal axis
    matrix = [[None for i in range(len(seq2)+1)] for j in range(len(seq1)+1)]

    # fill 0th row and 0th column with all 0
    for i in range(len(seq1)+1):
        matrix[i][0] = 0
    for i in range(len(seq2)+1):
        matrix[0][i] = 0

    print_2d_list(matrix)#test




#helper function to print matrix
def print_2d_list(matrix):
    for i in range(len(matrix)):
        print(matrix[i])



if __name__ == "__main__":
    seq1 = "AATCGTTCGAC"
    seq2 = "AACGTTTTCGGCA"
    smith_waterman(seq1, seq2)