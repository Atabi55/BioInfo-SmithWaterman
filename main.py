from smith_waterman import smith_waterman

def main():

    print("Smith-Waterman local alignment")
    print("Enter sequences to be aligned, leave blank for default sequence")
    seq1 = input("Sequence 1: ")
    if seq1 == "":
        seq1 = "AGTCCTGACGT"
    seq2 = input("Sequence 2: ")
    if seq2 == "":
        seq2 = "AGCTGACCGTGA"

    use_blosum_input = input("Do you want to use blosum62?(y/n): ")

    if use_blosum_input.lower() == 'y':
        gap = int(input("Enter gap penalty (commonly -8): "))
        smith_waterman(seq1, seq2, gap=gap, use_blosum62=True)
    else:
        match = int(input("Enter match score: "))
        mismatch = int(input("Enter mismatch penalty: "))
        gap = int(input("Enter gap penalty: "))

        smith_waterman(seq1, seq2, match=match, mismatch=mismatch, gap=gap, use_blosum62=False)

if __name__ == "__main__":
    main()
