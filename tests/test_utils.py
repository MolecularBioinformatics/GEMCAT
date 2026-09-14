from pathlib import Path

import cobra
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from fixtures import to_dense

from gemcat import utils

R_TOLERANCE = 10**-5
A_TOLERANCE = 10**-3


@pytest.fixture
def models():
    modelpath = Path("./tests/test_models/")
    modelpaths = {
        "mini": modelpath / "mini.xml",
        "mini_redox": modelpath / "mini_redox.xml",
        "mini_reversible": modelpath / "mini_reversible.xml",
    }
    models = {
        name: cobra.io.read_sbml_model(p.as_posix()) for name, p in modelpaths.items()
    }
    return models


def test_split_matrix_pos_neg():
    S = np.array(
        [
            [
                0,
                1,
                1,
                0,
                1,
                -2,
            ],
            [
                0,
                0,
                -1,
                2,
                -5,
                3,
            ],
            [
                -5,
                2,
                1,
                0,
                0,
                0,
            ],
        ]
    )
    result_pos, result_neg = utils.split_matrix_pos_neg(sp.csr_array(S.astype(float)))
    expected_pos = np.array(
        [
            [
                0,
                1,
                1,
                0,
                1,
                0,
            ],
            [
                0,
                0,
                0,
                2,
                0,
                3,
            ],
            [
                0,
                2,
                1,
                0,
                0,
                0,
            ],
        ]
    )
    expected_neg = np.array(
        [
            [
                0,
                0,
                0,
                0,
                0,
                -2,
            ],
            [
                0,
                0,
                -1,
                0,
                -5,
                0,
            ],
            [
                -5,
                0,
                0,
                0,
                0,
                0,
            ],
        ]
    )
    assert np.allclose(to_dense(result_pos), expected_pos, rtol=R_TOLERANCE)
    assert np.allclose(to_dense(result_neg), expected_neg, rtol=R_TOLERANCE)


def test_annotate():
    A = np.array([0, 1, 2])
    mets = ["m1", "m2", "m3"]
    result = utils.annotate_scores(A, mets)
    assert isinstance(result, pd.Series)
    assert len(result) == 3
    assert np.allclose(result.values, A, rtol=R_TOLERANCE)
    assert (result.index == mets).all()


def test_stoich_matrix_mini(models):
    model = models["mini"]
    result = utils.get_stoich_matrix_from_cobra(model)
    expected = np.array(
        [
            [
                -1,
                -1,
                0,
                0,
            ],
            [
                1,
                0,
                -1,
                0,
            ],
            [
                0,
                1,
                0,
                -1,
            ],
            [
                0,
                0,
                1,
                1,
            ],
        ]
    )
    assert np.allclose(to_dense(result), to_dense(expected), rtol=R_TOLERANCE)


def test_make_unidirectional():
    S = np.array(
        [
            [
                5,
                7,
                7,
                0,
                0,
                0,
                -3,
                0,
                -1,
            ],
            [
                0,
                0,
                0,
                2,
                0,
                -5,
                0,
                1,
                0,
            ],
            [
                0,
                -5,
                0,
                0,
                2,
                0,
                0,
                -8,
                0,
            ],
            [
                -4,
                0,
                0,
                -1,
                0,
                0,
                0,
                0,
                9,
            ],
            [
                0,
                -2,
                0,
                0,
                0,
                1,
                0,
                8,
                0,
            ],
        ]
    )
    expected = np.array(
        [
            [
                +5,
                +7,
                +7,
                00,
                00,
                00,
                -3,
                00,
                -1,
                00,
                00,
                +3,
                +1,
            ],
            [
                00,
                00,
                00,
                +2,
                00,
                -5,
                00,
                +1,
                00,
                -2,
                00,
                00,
                00,
            ],
            [
                00,
                -5,
                00,
                00,
                +2,
                00,
                00,
                -8,
                00,
                00,
                -2,
                00,
                00,
            ],
            [
                -4,
                00,
                00,
                -1,
                00,
                00,
                00,
                00,
                +9,
                +1,
                00,
                00,
                -9,
            ],
            [
                00,
                -2,
                00,
                00,
                00,
                +1,
                00,
                +8,
                00,
                00,
                00,
                00,
                00,
            ],
        ]
    )
    reversibilities = [0, 0, 0, 1, 1, 0, 1, 0, 1]
    reversibilities = [bool(r) for r in reversibilities]
    result = utils.make_unidirectional(sp.csc_array(S.astype(float)), reversibilities)
    assert np.allclose(to_dense(result), expected, rtol=R_TOLERANCE)


def test_make_unidirectional_wrong_type_int():
    S = np.array(
        [
            [
                5,
                7,
                7,
                0,
                0,
                0,
                -3,
                0,
                -1,
            ],
            [
                0,
                0,
                0,
                2,
                0,
                -5,
                0,
                1,
                0,
            ],
            [
                0,
                -5,
                0,
                0,
                2,
                0,
                0,
                -8,
                0,
            ],
            [
                -4,
                0,
                0,
                -1,
                0,
                0,
                0,
                0,
                9,
            ],
            [
                0,
                -2,
                0,
                0,
                0,
                1,
                0,
                8,
                0,
            ],
        ]
    )
    expected = np.array(
        [
            [
                +5,
                +7,
                +7,
                00,
                00,
                00,
                -3,
                00,
                -1,
                00,
                00,
                +3,
                +1,
            ],
            [
                00,
                00,
                00,
                +2,
                00,
                -5,
                00,
                +1,
                00,
                -2,
                00,
                00,
                00,
            ],
            [
                00,
                -5,
                00,
                00,
                +2,
                00,
                00,
                -8,
                00,
                00,
                -2,
                00,
                00,
            ],
            [
                -4,
                00,
                00,
                -1,
                00,
                00,
                00,
                00,
                +9,
                +1,
                00,
                00,
                -9,
            ],
            [
                00,
                -2,
                00,
                00,
                00,
                +1,
                00,
                +8,
                00,
                00,
                00,
                00,
                00,
            ],
        ]
    )
    reversibilities = [0, 0, 0, 1, 1, 0, 1, 0, 1]
    with pytest.raises(TypeError):
        result = utils.make_unidirectional(S, reversibilities)


def test_get_reversibilities_mini(models):
    model = models["mini"]
    expected = [False, False, False, False]
    result = utils.get_reversibilities(model)
    assert result == expected


def test_get_reversibilities_reversible(models):
    model = models["mini_reversible"]
    expected = [True, False, False, False]
    result = utils.get_reversibilities(model)
    assert result == expected


def test_get_metabolite_ids(models):
    model = models["mini"]
    result = utils.get_metabolite_ids(model)
    expected = ["A", "B", "C", "D"]
    assert result == expected


def test_geometric_mean():
    nums = [1, 2, 3, 4, 5, 6, 7]
    result = utils.geometric_mean(*nums)
    expected = 3.3800
    assert np.isclose(result, expected, rtol=R_TOLERANCE)


def test_geometric_mean_zero():
    nums = [0, 0, 0]
    result = utils.geometric_mean(*nums)
    expected = 0.0
    assert np.isclose(result, expected, rtol=R_TOLERANCE)


def test_geometric_mean_one():
    nums = [17]
    result = utils.geometric_mean(*nums)
    expected = 17.0
    assert np.isclose(result, expected, rtol=R_TOLERANCE)


def test_multiply_list():
    testcase = [1, 2, 3, 4, 5]
    result = utils.multiply(testcase)
    expected = 1 * 2 * 3 * 4 * 5
    print(result)
    print(expected)
    assert np.isclose(result, expected, rtol=R_TOLERANCE)


def test_multiply_single():
    testcase = 3
    result = utils.multiply(testcase)
    expected = 3
    assert np.isclose(result, expected, rtol=R_TOLERANCE)


def test_multiply_numbers():
    testcase = [1, 2, 3, 4, 5]
    result = utils.multiply(*testcase)
    expected = 1 * 2 * 3 * 4 * 5
    assert np.isclose(result, expected, rtol=R_TOLERANCE)
