"""Calculate the moment of inertia"""

import numpy as np


def principal_moments(masses_amu: np.ndarray, coords_A: np.ndarray) -> np.ndarray:
    """Principal moments of inertia I_A <= I_B <= I_C.

    Parameters
    -----------
    masses_amu : np.ndarray(N)
        Atomic masses in amu, N is the number of atoms

    coords_A : np.ndarray(N, 3)
        Cartesian coordinates of the molecule in A

    Returns
    ----------
    np.ndarray(I_A, I_B, I_C)
        Principal moments of inertia (from smallest to largest) in amu A^2

    References
    ----------
    D. A. McQuarrie, Statistical Mechanics (2000), ch. 5.
    """
    CM = np.average(coords_A, axis=0, weights=masses_amu)
    centered_coords_A = coords_A - CM
    mk_ri_rj = np.einsum(
        "i,ij,ik->jk", masses_amu, centered_coords_A, centered_coords_A
    )
    sum_m_r_2 = np.trace(mk_ri_rj) * np.eye(3)
    I = sum_m_r_2 - mk_ri_rj
    I_A_I_B_I_C = np.linalg.eigvalsh(I)

    return I_A_I_B_I_C
