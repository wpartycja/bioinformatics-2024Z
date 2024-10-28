from matrix import Matrix


class NeedlemanWunsch(Matrix):

    def initialize_matrix(self):
        """
        Initialize the scoring matrix with gap penalties.
        """

        super().initialize_matrix()

        # Initialize first column (gap in sequence1)
        for row in range(1, len(self.sequence1) + 1):
            self.matrix[row][0].value = self.gap_penalty * row
            self.matrix[row][0].max_from = "up"

        # Initialize first row (gap in sequence2)
        for column in range(1, len(self.sequence2) + 1):
            self.matrix[0][column].value = self.gap_penalty * column
            self.matrix[0][column].max_from = "left"

    def fill_matrix(self):
        """
        Fill the scoring matrix according to the Needleman-Wunsch algorithm.
        """

        for x in range(1, len(self.sequence1) + 1):
            for y in range(1, len(self.sequence2) + 1):

                score_up, score_diagonal, score_left = self.get_scores(x, y)

                # Select the maximum score and its direction
                max_score = max(score_diagonal, score_up, score_left)
                if max_score == score_diagonal:
                    self.matrix[x][y].max_from = "diagonal"
                elif max_score == score_up:
                    self.matrix[x][y].max_from = "up"
                else:
                    self.matrix[x][y].max_from = "left"

                self.matrix[x][y].value = max_score

    def traceback(self):
        """
        Perform traceback to get the optimal alignment.
        :return: Aligned sequences as a tuple.
        """
        x, y = len(self.sequence1), len(self.sequence2)

        return self.traceback_algortihm(x, y, "nw")

    def calculate_alignment_score(self) -> int:
        """
        Calculate the final alignment score based on the bottom-right cell of the matrix.
        :return: Final alignment score.
        """
        # Ensure the matrix is filled before calculating the score
        if self.matrix is None:
            self.perform_algorithm()  

        return self.matrix[len(self.sequence1)][len(self.sequence2)].value