from math import pi

import pytest

import qmthermo
from qmthermo import constants


def test_version():
    assert qmthermo.__version__ == "0.1.0"

def test_constants():
    T = 298.15
    assert constants.h / (2 * pi) == pytest.approx(1.054571817e-34, abs=1e-37)
    assert constants.kb * constants.Na == pytest.approx(8.314462618, abs=1e-3)
    assert constants.Eh_to_kcal_mol * 0.53108470 == pytest.approx(333.26098, abs=1e-3)
