import numpy as np
import pytest

from qmthermo.thermo import inertia

mass_vector = np.loadtxt("tests/mass_vector.txt", usecols=(2, 3, 4, 5))
masses_amu = mass_vector[:, 0]
coords_A = mass_vector[:, 1:4]


def test_read():
    assert sorted(
        inertia.principal_moments(masses_amu, coords_A), reverse=True
    ) == pytest.approx((2991.8376, 2283.13751, 2095.16738), abs=1e-3)


masses_H = np.array((1.007825032, 1.007825032))
coords_H = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 7.4]])


def test_physics():
    assert np.min(inertia.principal_moments(masses_H, coords_H)) == pytest.approx(
        0.0, abs=1e-5
    )
