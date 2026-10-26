"""Result type shared by all thermodynamic contribution"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Contribution:
    """Contributions to the thermodynamic functions for an ideal gas
    Enthalpy H in J/mol
    Entropy S in J/(mol K)
    Heat capacity Cv in J/(mol K)
    """

    H: float  # J mol^-1
    S: float  # J mol^-1 K^-1
    Cv: float  # J mol^-1 K^-1
