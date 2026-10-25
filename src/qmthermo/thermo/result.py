from dataclasses import dataclass


@dataclass
class Contribution:
    H: float   # J mol^-1
    S: float   # J mol^-1 K^-1
    Cv: float  # J mol^-1 K^-1
