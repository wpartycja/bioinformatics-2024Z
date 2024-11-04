from needlman_wunsch import NeedlemanWunsch
from smith_waterman import SmithWaterman


if __name__ == "__main__":
    # Needleman-Wunsch alignment
    print("----- Needleman Wunsch -----\n")

    # Example 1
    substitution_matrix_file = "./example_matrix.csv" 
    nw = NeedlemanWunsch("TTCA", "TGCTCG", 3, substitution_matrix_file)
    alignment = nw.perform_algorithm()
    nw.print_matrix()

    print("\nAligned Sequences:")
    for seq1, seq2, score in alignment:
        print(f"Score: {score}")
        print(f"Sequence 1: {seq1}")
        print(f"Sequence 2: {seq2}")
    
    print("\n\n")
    
    # Example 2
    nw = NeedlemanWunsch("TATA", "ATAT", 2, substitution_matrix_file)
    alignment = nw.perform_algorithm()
    nw.print_matrix()

    print("\nAligned Sequences:")
    for seq1, seq2, score in alignment:
        print(f"Score: {score}")
        print(f"Sequence 1: {seq1}")
        print(f"Sequence 2: {seq2}")
    
    print("\n\n")

    
    # Smith-Waterman alignment
    print("----- Smith-Waterman -----\n")

    # Example 1
    sw = SmithWaterman("TTCATA", "TGCTCGTA", 1, "./example_matrix.csv")
    alignment = sw.perform_algorithm()
    sw.print_matrix() 
    
    print(alignment)

    print("\nAligned Sequences:")
    for seq1, seq2, score in alignment:
        print(f"Score: {score}")
        print(f"Sequence 1: {seq1}")
        print(f"Sequence 2: {seq2}")

    print("\n\n")


    # Example 2
    sw = SmithWaterman("TATA", "ATAT", 1, "example_matrix.csv")
    alignment = sw.perform_algorithm()
    sw.print_matrix()  

    print("\nAligned Sequences:")
    for seq1, seq2, score in alignment:
        print(f"Score: {score}")
        print(f"Sequence 1: {seq1}")
        print(f"Sequence 2: {seq2}")

    print("\n\n")