from matrix import Matrix

class NeedlemanWunsch(Matrix):
    """
    For **global** alignment using the Needleman-Wunsch algorithm.
    """

    def initialize_matrix(self):
        """
        Initialize the scoring matrix with gap penalties.
        """
        super().initialize_matrix()

        # Initialize first column (gap in sequence1)
        for row in range(1, len(self.sequence1) + 1):
            self.matrix[row][0].value = self.gap_penalty * row
            self.matrix[row][0].max_from = ["up"]

        # Initialize first row (gap in sequence2)
        for column in range(1, len(self.sequence2) + 1):
            self.matrix[0][column].value = self.gap_penalty * column
            self.matrix[0][column].max_from = ["left"]

    def fill_matrix(self):
        """
        Fill the scoring matrix according to the Needleman-Wunsch algorithm.
        """
        for x in range(1, len(self.sequence1) + 1):
            for y in range(1, len(self.sequence2) + 1):
                score_up, score_diagonal, score_left = self.get_scores(x, y)

                # Select the maximum score and all directions that achieve it
                max_score = max(score_diagonal, score_up, score_left)
                directions = []
                if max_score == score_diagonal:
                    directions.append("diagonal")
                if max_score == score_up:
                    directions.append("up")
                if max_score == score_left:
                    directions.append("left")

                self.matrix[x][y].value = max_score
                self.matrix[x][y].max_from = directions

    def traceback(self):
        """
        Perform traceback to get all optimal alignments for global alignment.
        :return: List of tuples containing aligned sequences and their scores.
        """
        alignments = []
        x, y = len(self.sequence1), len(self.sequence2)
        self._traceback_recursive(x, y, [], [], alignments)
        
        curr_n = 0
        final_alignments = []
        for aligned_seq1, aligned_seq2 in alignments:
            if curr_n == self.n:
                break
            final_alignments.append((aligned_seq1, aligned_seq2, self.matrix[len(self.sequence1)][len(self.sequence2)].value))
            curr_n += 1
        
        return final_alignments

    def _traceback_recursive(self, x, y, aligned_seq1, aligned_seq2, alignments):
        """
        Helper method to recursively explore all paths in traceback.
        """
        # Base case: reach the top-left corner
        if x == 0 and y == 0:
            alignments.append((''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))))
            return
        
        # Recursive cases: explore all possible directions
        if "diagonal" in self.matrix[x][y].max_from:
            self._traceback_recursive(x - 1, y - 1, aligned_seq1 + [self.sequence1[x - 1]], aligned_seq2 + [self.sequence2[y - 1]], alignments)
        
        if "up" in self.matrix[x][y].max_from:
            self._traceback_recursive(x - 1, y, aligned_seq1 + [self.sequence1[x - 1]], aligned_seq2 + ['-'], alignments)
        
        if "left" in self.matrix[x][y].max_from:
            self._traceback_recursive(x, y - 1, aligned_seq1 + ['-'], aligned_seq2 + [self.sequence2[y - 1]], alignments)
