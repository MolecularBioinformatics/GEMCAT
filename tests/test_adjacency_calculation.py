from random import randint

import numpy as np
from fixtures import A_examples, S_examples

from gemcat import adjacency_transformation as at

R_TOLERANCE = 10**-5

reversibilities_linear = [False for i in range(6)]
expression_linear = np.ones((1, 6))

reversibilities_circular = [False for i in range(6)]
expression_circular = np.ones((1, 6))

reversibilities_bidir_linear = [False for i in range(6)]
expression_bidir_linear = np.ones((1, 6))

expression_rand_model = np.ones((1, 10))
reversibilities_rand_model = [False for i in range(10)]

reversibilities_complex = [False for i in range(18)]
expression_complex = np.ones((1, 18))


def test_AT_linear_ATHalf(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["linear"], reversibilities_linear, expression_linear, at.ATHalfStoich
    )
    assert np.allclose(result, A_examples["linear"], rtol=R_TOLERANCE)


def test_AT_linear_ATFull(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["linear"], reversibilities_linear, expression_linear, at.ATFullStoich
    )
    assert np.allclose(result, A_examples["linear"], rtol=R_TOLERANCE)


def test_AT_linear_ATPureAdj(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["linear"],
        reversibilities_linear,
        expression_linear,
        at.ATPureAdjacency,
    )
    assert np.allclose(result, A_examples["linear"], rtol=R_TOLERANCE)


def test_AT_circular_ATHalf(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["circular"],
        reversibilities_circular,
        expression_circular,
        at.ATHalfStoich,
    )
    assert np.allclose(result, A_examples["circular"], rtol=R_TOLERANCE)


def test_AT_circular_ATFull(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["circular"],
        reversibilities_circular,
        expression_circular,
        at.ATFullStoich,
    )
    assert np.allclose(result, A_examples["circular"], rtol=R_TOLERANCE)


def test_AT_circular_ATPureAdj(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["circular"],
        reversibilities_circular,
        expression_circular,
        at.ATPureAdjacency,
    )
    assert np.allclose(result, A_examples["circular"], rtol=R_TOLERANCE)


def test_AT_bidir_linear_ATHalf(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["bidir_linear"],
        reversibilities_bidir_linear,
        expression_bidir_linear,
        at.ATHalfStoich,
    )
    assert np.allclose(result, A_examples["bidir_linear"], rtol=R_TOLERANCE)


def test_AT_bidir_linear_ATFull(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["bidir_linear"],
        reversibilities_bidir_linear,
        expression_bidir_linear,
        at.ATFullStoich,
    )
    assert np.allclose(result, A_examples["bidir_linear"], rtol=R_TOLERANCE)


def test_AT_bidir_linear_ATPureAdj(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["bidir_linear"],
        reversibilities_bidir_linear,
        expression_bidir_linear,
        at.ATPureAdjacency,
    )
    assert np.allclose(result, A_examples["bidir_linear"], rtol=R_TOLERANCE)


def test_AT_complex_ATHalf(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["complex"],
        reversibilities_complex,
        expression_complex,
        at.ATHalfStoich,
    )
    assert np.allclose(A_examples["complex"], result, rtol=R_TOLERANCE)


def test_AT_complex_ATFull(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["complex"],
        reversibilities_complex,
        expression_complex,
        at.ATFullStoich,
    )
    assert np.allclose(A_examples["complex"], result, rtol=R_TOLERANCE)


def test_AT_complex_ATPureAdj(S_examples, A_examples):
    result = at.run_adjacencies_normalize(
        S_examples["complex"],
        reversibilities_complex,
        expression_complex,
        at.ATPureAdjacency,
    )
    assert np.allclose(A_examples["complex"], result, rtol=R_TOLERANCE)


def test_A_rows_zero_or_one_ATHalf():
    rand_vals = [randint(-5, 5) for i in range(100)]
    S = np.array(rand_vals).reshape(10, 10)
    A = at.run_adjacencies_normalize(
        S,
        reversibilities_rand_model,
        expression_rand_model,
        at.ATHalfStoich,
    )
    sums = A.sum(axis=1)
    for i in sums:
        assert np.isclose(i, 0.0, rtol=R_TOLERANCE) or np.isclose(
            i, 1.0, rtol=R_TOLERANCE
        )


def test_A_rows_zero_or_one_ATFull():
    rand_vals = [randint(-5, 5) for i in range(100)]
    S = np.array(rand_vals).reshape(10, 10)
    A = at.run_adjacencies_normalize(
        S,
        reversibilities_rand_model,
        expression_rand_model,
        at.ATFullStoich,
    )
    sums = A.sum(axis=1)
    for i in sums:
        assert np.isclose(i, 0.0, rtol=R_TOLERANCE) or np.isclose(
            i, 1.0, rtol=R_TOLERANCE
        )


def test_A_rows_zero_or_one():
    rand_vals = [randint(-5, 5) for i in range(100)]
    S = np.array(rand_vals).reshape(10, 10)
    A = at.run_adjacencies_normalize(
        S,
        reversibilities_rand_model,
        expression_rand_model,
        at.ATPureAdjacency,
    )
    sums = A.sum(axis=1)
    for i in sums:
        assert np.isclose(i, 0.0, rtol=R_TOLERANCE) or np.isclose(
            i, 1.0, rtol=R_TOLERANCE
        )
