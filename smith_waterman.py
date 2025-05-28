




def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2, use_blosum62=False):
    #Inputs: 2 DNA sequences, match, mismatch, gap scores, use_blosum62 boolean
    #Default values: match=2, mismatch=-1, gap=-2, use_blosum62=False
    #Output: List of all local alignments with inserted gaps


    #matrix: 2-dimensional list to hold score of each cell
    #seq1 in vertical axis, seq2 in horizontal axis
    #initialy filled with -1 since this value is not expected to occur in the matrix
    matrix = [[-1 for i in range(len(seq2)+1)] for j in range(len(seq1)+1)]

    # fill 0th row and 0th column with all 0
    for i in range(len(seq1)+1):
        matrix[i][0] = 0
    for i in range(len(seq2)+1):
        matrix[0][i] = 0

    print("Empty list:")
    print_2d_list(matrix)#test

    #arrows_matrix: 2-dimensional list, same size as matrix
    #holds the direction from which the calculated score comes
    #possible values: [], [diagonal], [up], [left], any combination of these directions in a list
    arrows_matrix = [[[] for i in range(len(seq2) + 1)] for j in range(len(seq1) + 1)]


    #loop through all cells in the matrix to calculate scores
    #also record from which direction the score comes
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            #recurrence: stores all 4 values of the recurrence relation for a particular cell
            recurrence = []
            #score from diagonal:
            recurrence.append(matrix[i][j] + match_mismatch_score(seq1[i], seq2[j], match, mismatch, use_blosum62))
            recurrence.append(matrix[i][j+1] + gap)#score from up (gap inserted to seq2)
            recurrence.append(matrix[i+1][j] + gap)#score from left (gap inserted to seq1)
            recurrence.append(0)

            #max of 4 recurrence values is the new value of this cell
            max_val = max(recurrence)
            matrix[i+1][j+1] = max_val

            #record the direction from which the max value comes
            if max_val != 0:
                if recurrence[0] == max_val:#diagonal
                    arrows_matrix[i+1][j+1].append("diagonal")
                if recurrence[1] == max_val:
                    arrows_matrix[i+1][j+1].append("up")
                if recurrence[2] == max_val:
                    arrows_matrix[i+1][j+1].append("left")

    print("All values of the matrix:")
    print_2d_list(matrix)  # test
    print("Directions matrix:")
    print_2d_list(arrows_matrix)  # test




#returns the score for comparing only 2 characters, blosum62 is used if specified.
def match_mismatch_score(c1, c2, match, mismatch, use_blosum62=False):
    if use_blosum62:
        #TODO implement later
        #pass
        return None
    else:
        if c1 == c2:
            return match
        else:
            return mismatch


#helper function to print matrix
def print_2d_list(x):
    for i in range(len(x)):
        print(x[i])



if __name__ == "__main__":
    seq1 = "AATCGTTCGAC"
    seq2 = "AACGTTTTCGGCA"
    smith_waterman(seq1, seq2)