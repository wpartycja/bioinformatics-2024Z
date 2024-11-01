from cell import Cell
from abc import ABC, abstractmethod
from tabulate import tabulate
import csv
import os

OUTPUT_PATH = './output'


class Matrix(ABC):
    """
    Virtual base class for alignment algorithms (Needleman-Wunsch, Smith-Waterman).
    """

    def __init__(
        self,
        sequence1: str,
        sequence2: str,
        n: int,
        filepath: str,
        gap_penalty: int = -2,
        output_filename: str = None
    ):
        """
        Initializes the base matrix for alignment algorithms.

        :param sequence1: the first sequence (column)
        :param sequence2: the second sequence (row)
        :param n: maximum number of optimal alignments
        :param filepath: filepath to the substitution matrix in CSV format
        :param gap_penalty: penalty for introducing gaps in the alignment (default -2)
        :param output_filename: name of the file to save the alignment output
        """
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        
        self.n = 1 if self.__class__.__name__ == 'SmithWaterman' else n
        self.gap_penalty = gap_penalty
        
        if not os.path.exists(OUTPUT_PATH):
            os.makedirs(OUTPUT_PATH)
        
        self.output_filename = output_filename if output_filename else f"{OUTPUT_PATH}/{self.__class__.__name__.lower()}_{self.sequence1}_{self.sequence2}.txt"
        self.substitution_matrix = self.load_matrix(filepath)

        self.matrix = None

    @staticmethod
    def load_matrix(path: str) -> dict:
        """
        Loads a substitution matrix that includes match, mismatch, and gap penalties from a CSV file.

        Parameters:
            - path (str): filepath to the substitution matrix in CSV format

        Returns:
            - matrix_ (dict): substitution matrix in a form of a nested dictionary where each key is a nucleotide,
              and each value is a dictionary mapping other nucleotides to their scores.
        """
        matrix_ = {}
        with open(path, 'r') as file:
            r = csv.reader(file)
            nucleotides1 = [header.strip() for header in next(r)[1:]]
            for row in r:
                nucleotide = row[0].strip()
                values = list(map(int, [x.strip() for x in row[1:]]))
                matrix_[nucleotide] = dict(zip(nucleotides1, values))
        return matrix_

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
        """
        Computes the scores for the cell at position (x, y) based on the substitution matrix.

        :param x: Row index in the scoring matrix (for sequence1)
        :param y: Column index in the scoring matrix (for sequence2)
        :return: A list containing the scores for up, diagonal, and left moves.
        """
        # Score for the diagonal move using the substitution matrix
        score_diagonal = (
            self.matrix[x - 1][y - 1].value +
            self.substitution_matrix[self.sequence1[x - 1]].get(self.sequence2[y - 1], self.gap_penalty)
        )

        # Score for the up move
        score_up = self.matrix[x - 1][y].value + self.gap_penalty

        # Score for the left move
        score_left = self.matrix[x][y - 1].value + self.gap_penalty

        return score_up, score_diagonal, score_left


    @abstractmethod
    def fill_matrix(self):
        """
        Fill the scoring matrix based on the specific algorithm (to be implemented by derived classes).
        """
        pass

    def traceback_algorithm(self, x, y, alg) -> list[str]:
        aligned_sequence1 = []
        aligned_sequence2 = []

        if alg == "nw":
            condition = lambda x, y: x > 0 or y > 0
        elif alg == "sw":
            condition = lambda x, y: x > 0 and y > 0 and self.matrix[x][y].value > 0
        else:
            raise ValueError("Invalid algorithm name for traceback!")

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
        alignments = self.traceback()
        self.save_output(alignments)
        
        return alignments
        

    def save_output(self, alignments):
        """
        Save the formatted output of the alignments to a file.
        """
        try:
            with open(self.output_filename, 'w') as f:
                for idx, (seq1, seq2, score) in enumerate(alignments):
                    f.write(f"Global alignment no. {idx + 1}:\n")
                    f.write(f"{seq1}\n")
                    f.write(f"{seq2}\n")
                    f.write(f"Score: {score}\n\n")
            print(f"Alignments saved successfully to {self.output_filename}.")
        except Exception as e:
            print(f"An error occurred while saving the alignments: {e}")
