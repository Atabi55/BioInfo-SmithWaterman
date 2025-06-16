import pandas as pd

#creating blosum62 matrix for future use
blosum62 = pd.read_csv( 'blosum62.dat', index_col=0, header=0, sep=r'\s+' )

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

    #print("Empty list:")
    #print_2d_list(matrix)#test

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
            
            

    #print("All values of the matrix:")
    #print_2d_list(matrix)  # test
    #print("Directions matrix:")
    #print_2d_list(arrows_matrix)  # test


    #Traceback starts at cell(s) with highest score
    #Find the highest score in matrix
    high_score = max([max(x) for x in matrix])

    alignment_results = []

    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] == high_score:
                traceback(seq1, seq2, arrows_matrix, matrix, i, j, "", "", alignment_results)
    
    print(f"Highest score: {high_score}")
    print(f"Number of optimal alignments: {len(alignment_results)}")
    for idx, (a1, a2) in enumerate(alignment_results, start=1):

        #making the "|" bars between two aligned subsequences
        bars = ["|" if a1[i] == a2[i] else " " for i in range(len(a1))]

        print(f"\nAlignment {idx}:")
        print(a1)
        print("".join(bars))
        print(a2)

    return alignment_results


    #count the number of highest value cells in each sublist (row) of matrix with matrix[i].count(high_score)
    #locate the index of each such cell in the list with matrix[i].index(high_score, start)
    #for each discovered cell with highest value, use the arrows_matrix to traceback until cell with value 0 reached
    #record a new pair of strings for each traceback, inserting sequence character or gap based on arrows_matrix
    #resulting strings are the final output.

def traceback(seq1, seq2, arrows_matrix, matrix, i, j, aligned1, aligned2, alignment_results):
    if matrix[i][j] == 0:
        reversed_al1 = ""
        reversed_al2 = ""
        for char in reversed(aligned1):
            reversed_al1 += char

        for char in reversed(aligned2):
            reversed_al2 += char

        alignment_results.append((reversed_al1, reversed_al2)) #reverse the strings and append to results
        return

    for direction in arrows_matrix[i][j]:
        if direction == "diagonal":
            traceback(seq1, seq2, arrows_matrix, matrix,
                      i-1, j-1,
                      aligned1 + seq1[i-1], aligned2 + seq2[j-1], alignment_results)
        elif direction == "up":
            traceback(seq1, seq2, arrows_matrix, matrix,
                      i-1, j,
                      aligned1 + seq1[i-1], aligned2 + "-", alignment_results)
        elif direction == "left":
            traceback(seq1, seq2, arrows_matrix, matrix,
                      i, j-1,
                      aligned1 + "-", aligned2 + seq2[j-1], alignment_results)




#returns the score for comparing only 2 characters, blosum62 is used if specified.
def match_mismatch_score(c1, c2, match, mismatch, use_blosum62=False):
    if use_blosum62:
        return blosum62[c1][c2]
    else:
        if c1 == c2:
            return match
        else:
            return mismatch




#helper function to print 2d list
def print_2d_list(x):
    for i in range(len(x)):
        print(x[i])



if __name__ == "__main__":
    #2 test cases
    seq1 = "AATCGTTCGAC"
    seq2 = "AACGTTTTCGGCA"
    smith_waterman(seq1, seq2)

    seq1 = "PHSWG"
    seq2 = "HGWAG"
    smith_waterman(seq1, seq2, gap=-8, use_blosum62=True)