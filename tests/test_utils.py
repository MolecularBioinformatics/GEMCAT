from pathlib import Path

import cobra
import numpy as np
import pandas as pd
import pytest

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


def test_l1_norm_1():
    a = np.array([1, 2, 3, 4, 5])
    l1 = utils._l1_norm(a)

    assert np.isclose(l1, 15.0, rtol=R_TOLERANCE)


def test_l1_norm_long():
    a = np.ones(100)
    l1 = utils._l1_norm(a)

    assert np.isclose(l1, 100.0, rtol=R_TOLERANCE)


def test_l1_norm_zeroes():
    a = np.zeros(100)
    l1 = utils._l1_norm(a)

    assert np.isclose(l1, 0.0, rtol=R_TOLERANCE)


def test_l1_empty():
    a = np.array([])
    with pytest.raises(ValueError):
        utils._l1_norm(a)


def test_get_ids():
    g1 = cobra.Gene("g1")
    g2 = cobra.Gene("g2")
    result = utils._get_ids([g1, g2])
    expected = ["g1", "g2"]

    assert result == expected


def test_get_n_reactions():
    S = np.array(
        [
            [
                0,
                1,
                1,
                0,
                1,
                2,
            ],
            [
                0,
                0,
                1,
                2,
                5,
                3,
            ],
            [
                5,
                2,
                1,
                0,
                0,
                0,
            ],
        ]
    )
    result = utils._get_n_reactions(S)
    expected = np.array(
        [
            4,
            4,
            3,
        ]
    )
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_get_total_stoich():
    S = np.array(
        [
            [
                0,
                1,
                1,
                0,
                1,
                2,
            ],
            [
                0,
                0,
                1,
                2,
                5,
                3,
            ],
            [
                5,
                2,
                1,
                0,
                0,
                0,
            ],
        ]
    )
    result = utils._get_total_stoich(S)
    expected = np.array(
        [
            5,
            11,
            8,
        ]
    )
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


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
    result_pos, result_neg = utils.split_matrix_pos_neg(S)
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
    assert np.allclose(result_pos, expected_pos, rtol=R_TOLERANCE)
    assert np.allclose(result_neg, expected_neg, rtol=R_TOLERANCE)


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
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_linearization(models):
    model = models["mini_reversible"]
    expected = np.array(
        [
            [-1, -1, 00, 00, +1],
            [+1, 00, -1, 00, -1],
            [00, +1, 00, -1, 00],
            [00, 00, +1, +1, 00],
        ]
    )
    result = utils._get_unidirectional_matrix(model)
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_linearization_already_linear(models):
    model = models["mini"]
    expected = utils.get_stoich_matrix_from_cobra(model)
    result = utils._get_unidirectional_matrix(model)
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_replace_zeroes_nan():
    A = np.array(
        [
            [
                0,
                0,
                np.nan,
            ],
            [
                0,
                np.nan,
                -1,
            ],
            [
                2,
                0,
                np.nan,
            ],
        ]
    )
    expected = np.array(
        [
            [
                0,
                0,
                0,
            ],
            [
                0,
                0,
                -1,
            ],
            [
                2,
                0,
                0,
            ],
        ]
    )
    result = utils._replace_zeroes(A)
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_replace_zeroes_inf():
    A = np.array(
        [
            [
                0,
                0,
                np.inf,
            ],
            [
                0,
                np.inf,
                -1,
            ],
            [
                -2,
                0,
                5,
            ],
        ]
    )
    expected = np.array(
        [
            [
                0,
                0,
                0,
            ],
            [
                0,
                0,
                -1,
            ],
            [
                -2,
                0,
                5,
            ],
        ]
    )
    result = utils._replace_zeroes(A)
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_replace_zeroes_neg_inf():
    A = np.array(
        [
            [
                0,
                0,
                -np.inf,
            ],
            [
                -np.inf,
                -10,
                -1,
            ],
            [
                -2,
                0,
                5,
            ],
        ]
    )
    expected = np.array(
        [
            [
                0,
                0,
                0,
            ],
            [
                0,
                -10,
                -1,
            ],
            [
                -2,
                0,
                5,
            ],
        ]
    )
    result = utils._replace_zeroes(A)
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_replace_zeroes_mixed():
    A = np.array(
        [
            [
                0,
                0,
                -np.inf,
            ],
            [
                np.inf,
                -10,
                -1,
            ],
            [
                -2,
                np.nan,
                5,
            ],
        ]
    )
    expected = np.array(
        [
            [
                0,
                0,
                0,
            ],
            [
                0,
                -10,
                -1,
            ],
            [
                -2,
                0,
                5,
            ],
        ]
    )
    result = utils._replace_zeroes(A)
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_calc_zscore():
    scores = pd.Series(
        {
            "A": 1,
            "B": 2,
            "C": 3,
            "D": 4,
            "E": 5,
        }
    )
    expected = np.array([-1.26491106, -0.63245553, 0.0, 0.63245553, 1.26491106])
    result = utils._calc_zscore(scores).values
    assert np.allclose(result, expected, atol=0.01)


def test_scale():
    scores = pd.Series(
        {
            "A": 1,
            "B": 2,
            "C": 3,
            "D": 4,
            "E": 5,
        }
    )
    result = utils._scale(scores)
    expected = pd.Series(
        {
            "A": 1 / 15,
            "B": 2 / 15,
            "C": 3 / 15,
            "D": 4 / 15,
            "E": 5 / 15,
        }
    )
    assert np.allclose(result, expected, atol=0.01)


def test_find_indeces():
    rxns = [
        "EX_ATP_123",
        "OF_NAD_abc",
        "R_12345",
        "R_EXTREME_BLA",
        "OFTEN_123",
        "EX_test",
    ]
    result = utils._find_indeces(rxns)
    expected = [2, 3, 4]
    assert result == expected


def test_is_exchange():
    rxns = ["EX_O2", "transport_o2", "EXTRACELLULAR", "OF_ATP", "OFTEN"]
    exchanges = [0, 3]
    for i, r in enumerate(rxns):
        result = utils._is_exchange(r)
        if i in exchanges:
            assert result
        else:
            assert not result


def test_get_subset_cols():
    indeces = [0, 1, 2, 8]
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
                5,
                7,
                7,
                -1,
            ],
            [
                0,
                0,
                0,
                0,
            ],
            [
                0,
                -5,
                0,
                0,
            ],
            [
                -4,
                0,
                0,
                9,
            ],
            [
                0,
                -2,
                0,
                0,
            ],
        ]
    )
    result = utils._get_subset_cols(S, indeces)
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_remove_exchanges():
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
    rxns = [
        "EX_O2",
        "OF_ATP",
        "EX_CO2",
        "transport_o2",
        "EXTRACELLULAR",
        "EXo2",
        "OFTEN",
        "E_lac",
        "EX_lac",
    ]
    expected = np.array(
        [
            [
                0,
                0,
                0,
                -3,
                0,
            ],
            [
                2,
                0,
                -5,
                0,
                1,
            ],
            [
                0,
                2,
                0,
                0,
                -8,
            ],
            [
                -1,
                0,
                0,
                0,
                0,
            ],
            [
                0,
                0,
                1,
                0,
                8,
            ],
        ]
    )
    result = utils._remove_exchanges(S, rxns)
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


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
    result = utils.make_unidirectional(S, reversibilities)
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


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


def test_get_reaction_ids(models):
    model = models["mini"]
    result = utils._get_reaction_ids(model)
    expected = ["R1", "R2", "R3", "R4"]
    assert result == expected


def test_get_metabolite_ids(models):
    model = models["mini"]
    result = utils.get_metabolite_ids(model)
    expected = ["A", "B", "C", "D"]
    assert result == expected


def test_make_row_vector():
    arr = np.array([1, 2, 3])
    result = utils.make_row_vector(arr)
    expected = np.array([[1, 2, 3]])
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_make_column_vector():
    arr = np.array([1, 2, 3])
    result = utils.make_column_vector(arr)
    expected = np.array(
        [
            [1],
            [2],
            [3],
        ]
    )
    assert np.allclose(result, expected, rtol=R_TOLERANCE)


def test_is_array():
    S = [
        [-1, 0, 0, 0, 0, 0],
        [1, -1, 0, 0, 0, 0],
        [0, 1, -1, 0, 0, 0],
        [0, 0, 1, -1, 0, 0],
        [0, 0, 0, 1, -1, 0],
        [0, 0, 0, 0, 1, 0],
    ]
    with pytest.raises(TypeError):
        utils._is_np_array(S)


def test_check_array_shape():
    Si = np.array(
        [
            [-1, 0, 0, 0, 0, 0],
            [1, -1, 0, 0, 0, 0],
            [0, 1, -1, 0, 0, 0],
            [0, 0, 1, -1, 0, 0],
            [0, 0, 0, 1, -1, 0],
            [0, 0, 0, 0, 1, 0],
        ]
    )
    Sj = np.array(
        [
            [1, -1, 0, 0, 0, 0],
            [0, 1, -1, 0, 0, 0],
            [0, 0, 1, -1, 0, 0],
            [0, 0, 0, 1, -1, 0],
            [0, 0, 0, 0, 1, 0],
        ]
    )
    with pytest.raises(ValueError):
        utils._check_array_shape(Si, Sj)


def test_is_all_ones_true():
    S = np.array(
        [
            [-1, 0, 0, 0, 0, 0],
            [1, -1, 0, 0, 0, 0],
            [0, 1, -1, 0, 0, 0],
            [0, 0, 1, -1, 0, 0],
            [0, 0, 0, 1, -1, 0],
            [0, 0, 0, 0, 1, 0],
        ]
    )
    result = utils._is_all_ones(S)
    expected = False
    assert result == expected


def test_is_all_ones_2d_true():
    S = np.ones((6, 6))
    result = utils._is_all_ones(S)
    expected = True
    assert result == expected


def test_is_all_ones_1d_true():
    S = np.ones((6))
    result = utils._is_all_ones(S)
    expected = True
    assert result == expected


def test_is_all_ones_1_5d_true():
    S = np.ones((6, 1))
    result = utils._is_all_ones(S)
    expected = True
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
