from math import e, log, pi

from qmthermo import constants
from qmthermo.thermo.result import Contribution


def translation(mass_amu: float, T: float, P_atm: float = 1.0) -> Contribution:
    """Translational contributions to enthalpy, entropy and heat capacity
    
    Parameters:
    -----------
    mass_amu : float
        Molecular mass im amu
        
    T : float
        Temperature in K
    
    P_atm : float
        Pressure in atm, 1.0 atm by default

    Returns:
    --------
    Contribution(H = ..., S = ..., Cv = ...)

    References
    ----------
    D. A. McQuarrie, Statistical Mechanics (2000), ch. 5.
    """
    m = mass_amu * constants.amu_in_kg
    l = ((2 * pi * m * constants.kb * T) / constants.h ** 2) ** (1/2)
    # #Helmholtz free energy
    # A_trans = -constants.R * T * log(q_trans)
    # Inner energy
    E_trans = 3/2 * constants.R * T
    # Entropy
    S_trans = constants.R * log(l ** 3 * constants.R * T * e ** (5/2) / (constants.Na * P_atm * constants.atm_to_Pa))
    # Enthalpy
    H_trans = E_trans + constants.R * T
    # Heat capacity
    Cv_trans = 3/2 * constants.R

    return Contribution(H=H_trans, S=S_trans, Cv=Cv_trans)
    