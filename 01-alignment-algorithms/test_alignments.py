import pytest
from needlman_wunsch import NeedlemanWunsch
from smith_waterman import SmithWaterman

# Fixtures for consistent parameters across tests
@pytest.fixture
def nw_params():
    return {
        "sequence1": "TTCATA",
        "sequence2": "TGCTCGTA",
        "gap_penalty": -6,
        "match_score": 5,
        "mismatch_penalty": -2
    }

@pytest.fixture
def sw_params():
    return {
        "sequence1": "TTCATA",
        "sequence2": "TGCTCGTA",
        "gap_penalty": -6,
        "match_score": 5,
        "mismatch_penalty": -2
    }

def test_needleman_wunsch_basic_alignment(nw_params):
    nw = NeedlemanWunsch(**nw_params)
    alignment = nw.perform_algorithm()
    
    expected_alignment1 = "T--TCATA"
    expected_alignment2 = "TGCTCGTA"
    
    assert alignment[0] == expected_alignment1, f"Expected {expected_alignment1}, got {alignment[0]}"
    assert alignment[1] == expected_alignment2, f"Expected {expected_alignment2}, got {alignment[1]}"

def test_smith_waterman_basic_alignment(sw_params):
    sw = SmithWaterman(**sw_params)
    alignment = sw.perform_algorithm()
    
    expected_alignment1 = "TCATA"
    expected_alignment2 = "TCGTA"
    
    assert alignment[0] == expected_alignment1, f"Expected {expected_alignment1}, got {alignment[0]}"
    assert alignment[1] == expected_alignment2, f"Expected {expected_alignment2}, got {alignment[1]}"

# Edge case tests for empty sequences
@pytest.mark.parametrize("sequence1, sequence2, expected_alignment1, expected_alignment2", [
    ("", "TGCTCGTA", "--------", "TGCTCGTA"),  
    ("TTCATA", "", "TTCATA", "------"),
    ("", "", "", "")
])
def test_needleman_wunsch_empty_sequences(sequence1, sequence2, expected_alignment1, expected_alignment2, nw_params):
    nw_params.update({"sequence1": sequence1, "sequence2": sequence2})
    nw = NeedlemanWunsch(**nw_params)
    alignment = nw.perform_algorithm()

    assert alignment[0] == expected_alignment1, f"Expected {expected_alignment1}, got {alignment[0]}"
    assert alignment[1] == expected_alignment2, f"Expected {expected_alignment2}, got {alignment[1]}"

# Edge case tests for empty sequences
@pytest.mark.parametrize("sequence1, sequence2, expected_alignment1, expected_alignment2", [
    ("", "TGCTCGTA", "", ""),  
    ("TTCATA", "", "", ""),    
    ("", "", "", "")
])
def test_smith_waterman_empty_sequences(sequence1, sequence2, expected_alignment1, expected_alignment2, sw_params):
    sw_params.update({"sequence1": sequence1, "sequence2": sequence2})
    sw = SmithWaterman(**sw_params)
    alignment = sw.perform_algorithm()

    assert alignment[0] == expected_alignment1, f"Expected {expected_alignment1}, got {alignment[0]}"
    assert alignment[1] == expected_alignment2, f"Expected {expected_alignment2}, got {alignment[1]}"

def test_needleman_wunsch_scoring(nw_params):
    nw = NeedlemanWunsch(**nw_params)
    nw.perform_algorithm()
    
    expected_score = 11 
    assert nw.calculate_alignment_score() == expected_score, f"Expected score {expected_score}, got {nw.calculate_alignment_score()}"

def test_smith_waterman_scoring(sw_params):
    sw = SmithWaterman(**sw_params)
    sw.perform_algorithm()
    
    expected_score = 18  # Replace with actual expected score for Smith-Waterman
    assert sw.calculate_alignment_score() == expected_score, f"Expected score {expected_score}, got {sw.calculate_alignment_score()}"

def test_needleman_wunsch_identical_sequences(nw_params):
    nw_params.update({"sequence1": "TTCATA", "sequence2": "TTCATA"})
    nw = NeedlemanWunsch(**nw_params)
    alignment = nw.perform_algorithm()

    assert alignment[0] == "TTCATA"
    assert alignment[1] == "TTCATA"

def test_smith_waterman_identical_sequences(sw_params):
    sw_params.update({"sequence1": "TTCATA", "sequence2": "TTCATA"})
    sw = SmithWaterman(**sw_params)
    alignment = sw.perform_algorithm()

    assert alignment[0] == "TTCATA"
    assert alignment[1] == "TTCATA"
