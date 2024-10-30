from needlman_wunsch import NeedlemanWunsch
from smith_waterman import SmithWaterman

if __name__ == "__main__":
    # Needlaman Wunsch alignment
    print("----- Needleman Wunsch -----\n")

    nw = NeedlemanWunsch("TTCATA", "TGCTCGTA", -6, 5, -2)
    alignment = nw.perform_algorithm()
    nw.print_matrix()

    print("\nAligned Sequences:")
    print(alignment)

    print("\n\n")
    
    nw = NeedlemanWunsch("TATA", "ATAT", -1, 1, -1)
    alignment = nw.perform_algorithm()
    nw.print_matrix()

    print("\nAligned Sequences:")
    print(alignment)
    

    print("\n\n")
    
    # Smith-Waterman alignment
    print("\n\n----- Smith-Waterman -----\n")

    sw = SmithWaterman("TTCATA", "TGCTCGTA", -6, 5, -2)
    alignment_sw = sw.perform_algorithm()
    sw.print_matrix()

    print("\nAligned Sequences:")
    print(alignment_sw)
    
    print("\n\n")
    
    sw = SmithWaterman("TATA", "ATAT", -1, 1, -1)
    alignment = sw.perform_algorithm()
    sw.print_matrix()

    print("\nAligned Sequences:")
    print(alignment)

