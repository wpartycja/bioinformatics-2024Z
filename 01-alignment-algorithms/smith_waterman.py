from matrix import Matrix


class SmithWaterman(Matrix):
    """
    Class that implements the Smith-Waterman algorithm (local alignment).
    """

    def fill_matrix(self):
        """
        Fill the matrix based on Smith-Waterman scoring (local alignment).
        """
        for x in range(1, len(self.sequence1) + 1):
            for y in range(1, len(self.sequence2) + 1):

                score_up, score_diagonal, score_left = self.get_scores(x, y)

                max_score = max(
                    0, score_diagonal, score_up, score_left
                )  # Smith-Waterman allows scores to be zero
                if max_score == score_diagonal:
                    self.matrix[x][y].max_from = "diagonal"
                elif max_score == score_up:
                    self.matrix[x][y].max_from = "up"
                elif max_score == score_left:
                    self.matrix[x][y].max_from = "left"
                else:
                    self.matrix[x][y].max_from = None  # No traceback if score is 0

                self.matrix[x][y].value = max_score

    def traceback(self) -> list[str]:
        """
        Perform the traceback to get the optimal alignment for Smith-Waterman (local alignment).
        """

        # Find the cell with the highest score for local alignment
        max_x, max_y = 0, 0
        max_value = 0
        for x in range(1, len(self.sequence1) + 1):
            for y in range(1, len(self.sequence2) + 1):
                if self.matrix[x][y].value > max_value:
                    max_value = self.matrix[x][y].value
                    max_x, max_y = x, y

        x, y = max_x, max_y

        return self.traceback_algortihm(x, y, "sw")


    def calculate_alignment_score(self) -> int:
        """
        Calculate the highest alignment score in the matrix for local alignment.
        
        :return: The highest score found in the matrix, representing the optimal local alignment score.
        """
        if self.matrix is None:
            self.perform_algorithm()  

        # Find the maximum value in the matrix
        max_value = 0
        for row in self.matrix:
            for cell in row:
                if cell.value > max_value:
                    max_value = cell.value

        return max_value