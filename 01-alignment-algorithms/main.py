from needlman_wunsch import NeedlemanWunsch
from smith_waterman import SmithWaterman

if __name__ == "__main__":
    # Needlaman Wunsch alignment
    print("----- Needleman Wunsch -----\n")

    nw = NeedlemanWunsch("TTCATA", "TGCTCGTA", -6, 5, -2)
    alignment = nw.perform_algorithm()
    nw.print_matrix()

    print("\nAligned Sequences:")
    print(alignment[0])
    print(alignment[1])
    
    nw_score = nw.calculate_alignment_score()
    
    print(f"Caluclated score: {nw_score}")

    # Smith-Waterman alignment
    print("\n\n----- Smith-Waterman -----\n")

    sw = SmithWaterman("TTCATA", "TGCTCGTA", -6, 5, -2)
    alignment_sw = sw.perform_algorithm()
    sw.print_matrix()

    print("\nAligned Sequences:")
    print(alignment_sw[0])
    print(alignment_sw[1])
    
    sw_score = sw.calculate_alignment_score()
    
    print(f"Caluclated score: {sw_score}")
