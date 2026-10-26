from math import log

import pytest

from qmthermo import constants
from qmthermo.thermo.translation import translation

T = 298.15
mass = 340.20931

def test_reproduce_H():
    assert translation(mass, T).H / constants.kcal_to_J == pytest.approx(1.48121, abs=1e-3)


def test_reproduce_S():
    assert translation(mass, T).S * 1000 / constants.kcal_to_J  == pytest.approx(43.36774, abs=1e-3)


def test_reproduce_Cv():
    assert translation(mass, T).Cv * 1000 / constants.kcal_to_J == pytest.approx(2.98081, abs=1e-3)    


def test_physics_S():
    lower_S = 43.36774 - constants.R * log(2) * 1000 / constants.kcal_to_J
    assert translation(mass, T, P_atm=2).S * 1000 / constants.kcal_to_J == pytest.approx(lower_S, abs=1e-3)


def test_physics_H():
    assert translation(mass, T, 2).H / constants.kcal_to_J == pytest.approx(1.48121, abs=1e-3)
