class Cell:
    """
    Represents a single cell in the Needleman-Wunsch matrix.

    :param value: The score of this cell.
    :param max_from: The direction from which the score was derived ('diagonal', 'up', 'left').
    """

    def __init__(self, value: int = 0, max_from: str = None):
        self.value = value
        self.max_from = max_from  # 'up', 'left', 'diagonal'
