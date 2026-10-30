"""Calculate the moment of inertia"""

import numpy as np


def principal_moments(masses_amu: np.ndarray, coords_A: np.ndarray) -> np.ndarray:
    """Principal moments of inertia (Ixx, Iyy, Izz).

    Parameters:
    -----------
    masses_amu : np.ndarray(N)
        Molecular masses in amu, N is the number of atoms

    coords_A : np.ndarray(N, 3)
        Cartesian coordinates of the molecule in A

    Returns:
    ----------
    np.ndarray(Ixx, Iyy, Izz)
        Principal moments of inertia (from smallest to largest) im amu A^2

    References
    ----------
    D. A. McQuarrie, Statistical Mechanics (2000), ch. 5.
    """
    CM = np.average(coords_A, axis=0, weights=masses_amu)
    centered_coords_A = coords_A - CM
    # I = np.zeros((3,3))
    # xi_xj = np.asarray([np.outer(coords_A_k, coords_A_k) for coords_A_k in centered_coords_A])
    mk_ri_rj = np.einsum(
        "i,ij,ik->jk", masses_amu, centered_coords_A, centered_coords_A
    )
    r_2 = np.asarray(
        [
            np.einsum("ij,ij->i", centered_coords_A, centered_coords_A)[i] * np.eye(3)
            for i in range(len(coords_A))
        ]
    )
    sum_m_r_2 = np.einsum("i,ijk->jk", masses_amu, r_2)
    I = sum_m_r_2 - mk_ri_rj
    Ixx_Iyy_Izz = np.linalg.eigvalsh(I)

    return Ixx_Iyy_Izz
