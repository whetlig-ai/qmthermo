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
    S_1 = translation(mass, T, 1).S 
    S_2 = translation(mass, T, 2).S 
    assert (S_1 - S_2) == pytest.approx(constants.R * log(2), abs=1e-3, rel=1e-12)


def test_physics_H():
    assert translation(mass, T, 2).H / constants.kcal_to_J == pytest.approx(1.48121, abs=1e-3)
