"""Translational contributions to
Enthalpy, Entropy and Heat capacity of theideal gas
"""

from math import log, pi, sqrt

from qmthermo import constants
from qmthermo.thermo.result import Contribution


def translation(mass_amu: float, T: float, P_atm: float = 1.0) -> Contribution:
    """Translational contributions to enthalpy, entropy and heat capacity

    Parameters
    -----------
    mass_amu : float
        Molecular mass in amu

    T : float
        Temperature in K

    P_atm : float
        Pressure in atm, 1.0 atm by default

    Returns
    --------
    Contribution(
    H = ... in J/mol,
    S = ... in J/(mol K) ,
    Cv = ... in J/(mol K))

    References
    ----------
    D. A. McQuarrie, Statistical Mechanics (2000), ch. 5.
    """
    m = mass_amu * constants.amu_in_kg
    Lambda = constants.h / sqrt(2 * pi * m * constants.kb * T)
    q = constants.kb * T / (P_atm * constants.atm_to_Pa * Lambda**3)
    # Inner energy
    E_trans = 3 / 2 * constants.R * T
    # Entropy
    S_trans = constants.R * (log(q) + 5 / 2)
    # Enthalpy
    H_trans = E_trans + constants.R * T
    # Heat capacity
    Cv_trans = 3 / 2 * constants.R

    return Contribution(H=H_trans, S=S_trans, Cv=Cv_trans)
