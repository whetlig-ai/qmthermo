from math import pi

import pytest

import qmthermo
from qmthermo import constants


def test_version():
    assert qmthermo.__version__ == "0.1.0"

def test_h_red():
    assert constants.h_red == pytest.approx(1.054571817e-34, rel=1e-9)

def test_R():
    T = 298.15
    assert constants.R == pytest.approx(8.314462618, abs=1e-3)
    assert constants.R * T / constants.kcal_to_J == pytest.approx(0.592485, abs=1e-3)
