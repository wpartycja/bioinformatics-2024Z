from cell import Cell
from abc import ABC, abstractmethod
from tabulate import tabulate


class Matrix(ABC):
    """
    Virtual base class for alignment algorithms (Needleman-Wunsch, Smith-Waterman).
    """

    def __init__(
        self,
        sequence1: str,
        sequence2: str,
        gap_penalty: int,
        match_score: int,
        mismatch_penalty: int,
    ):
        """
        Initializes the base matrix for alignment algorithms.

        :param sequence1: the first sequence (column)
        :param sequence2: the second sequence (row)
        :param gap_penalty: penalty for introducing gaps in the alignment
        :param match_score: score for matching characters
        :param mismatch_penalty: penalty for mismatching characters
        """
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.gap_penalty = gap_penalty
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty

        self.matrix = None

    def initialize_matrix(self):
        """
        Initialize the matrix with proper values (to be implemented by derived classes).
        """

        # Create a matrix with an extra row and column for the initial gaps
        self.matrix = [
            [Cell() for _ in range(len(self.sequence2) + 1)]
            for _ in range(len(self.sequence1) + 1)
        ]

    def get_scores(self, x: int, y: int) -> list[int]:
        if self.sequence1[x - 1] == self.sequence2[y - 1]:
            score_diagonal = self.matrix[x - 1][y - 1].value + self.match_score
        else:
            score_diagonal = self.matrix[x - 1][y - 1].value + self.mismatch_penalty

        score_up = self.matrix[x - 1][y].value + self.gap_penalty
        score_left = self.matrix[x][y - 1].value + self.gap_penalty

        return score_up, score_diagonal, score_left

    @abstractmethod
    def fill_matrix(self):
        """
        Fill the scoring matrix based on the specific algorithm (to be implemented by derived classes).
        """
        pass

    def traceback_algortihm(self, x, y, alg) -> list[str]:
        aligned_sequence1 = []
        aligned_sequence2 = []

        if alg == "nw":
            condition = lambda x, y: x > 0 or y > 0
        elif alg == "sw":
            condition = lambda x, y: x > 0 and y > 0 and self.matrix[x][y].value > 0
        else:
            ValueError("Invalid algorithm name for traceback!")

        while condition(x, y):
            if self.matrix[x][y].max_from == "diagonal":
                aligned_sequence1.append(self.sequence1[x - 1])
                aligned_sequence2.append(self.sequence2[y - 1])
                x -= 1
                y -= 1
            elif self.matrix[x][y].max_from == "up":
                aligned_sequence1.append(self.sequence1[x - 1])
                aligned_sequence2.append("-")
                x -= 1
            elif self.matrix[x][y].max_from == "left":
                aligned_sequence1.append("-")
                aligned_sequence2.append(self.sequence2[y - 1])
                y -= 1

        # Reverse the sequences as the traceback starts from the bottom-right
        aligned_sequence1.reverse()
        aligned_sequence2.reverse()

        return "".join(aligned_sequence1), "".join(aligned_sequence2)

    @abstractmethod
    def traceback(self):
        """
        Perform traceback to get the optimal alignment (to be implemented by derived classes).
        """
        pass

    def print_matrix(self):
        """
        Print the scoring matrix using tabulate for a nicer output.
        """
        headers = [" "] + list("-" + self.sequence2)

        table = []
        for x in range(len(self.sequence1) + 1):
            if x == 0:
                row_label = "-"
            else:
                row_label = self.sequence1[x - 1]

            row = [row_label] + [
                self.matrix[x][y].value for y in range(len(self.sequence2) + 1)
            ]
            table.append(row)

        print(tabulate(table, headers, tablefmt="fancy_grid"))

    def perform_algorithm(self) -> list[str]:
        """
        Performs the complete alignment process and returns the aligned sequences.
        """
        self.initialize_matrix()
        self.fill_matrix()
        return self.traceback()
