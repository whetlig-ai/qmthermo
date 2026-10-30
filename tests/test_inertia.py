from pathlib import Path

import numpy as np
import pytest

from qmthermo.thermo import inertia

DATA = Path(__file__).parent / "tiobu4_geometry.txt"

mass_vector = np.loadtxt(DATA, usecols=(2, 3, 4, 5))
masses_amu = mass_vector[:, 0]
coords_A = mass_vector[:, 1:4]


def test_tiobu4_principal_moments():
    assert sorted(
        inertia.principal_moments(masses_amu, coords_A), reverse=True
    ) == pytest.approx((2991.8376, 2283.13751, 2095.16738), abs=1e-4)


masses_H = np.array((1.007825032, 1.007825032))
coords_H = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
mu = masses_H[0] * masses_H[1] / (np.sum(masses_H))
expected = mu * 0.74**2


def test_H2_moments():
    assert sorted(
        inertia.principal_moments(masses_H, coords_H), reverse=True
    ) == pytest.approx((expected, expected, 0.0), abs=1e-5)
