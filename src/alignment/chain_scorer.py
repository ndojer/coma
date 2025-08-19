# from typing import List

# from src.alignment.alignment_position_scorer import AlignmentPositionScorer


class ChainScorer:
    def __init__(self, 
                 indelOpenPenalty: int,
                 indelExtPenalty: float,
                 minAlignmentScore: int):
        self.indelOpenPenalty = indelOpenPenalty
        self.indelExtPenalty = indelExtPenalty
        self.minAlignmentScore = minAlignmentScore
        self.ultimatePenalty = -10^6
