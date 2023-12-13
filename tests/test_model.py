import numpy as np
from fixtures import S_examples

from gemcat.model import Model

R_TOLERANCE = 10**-2  # why is the accuracy lower here than in test_Pagerank?


def test_basic_linear_model():
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
    rev = [False for i in range(6)]
    expected = [0.061, 0.112, 0.156, 0.193, 0.225, 0.252]
    met_names = ["M" + str(x) for x in range(6)]
    m = Model(S, met_names, rev)
    scores = m.calculate()
    assert np.allclose(m.scores, expected, rtol=R_TOLERANCE)


def test_circular_model_w_expression(S_examples):
    S = S_examples["circular"]
    rev = [False for i in range(6)]
    expected = np.ones(6) * 1 / 6
    met_names = ["M" + str(x) for x in range(6)]
    m = Model(S, met_names, rev)
    m.expression_vector = np.ones((1, 6)) * 2
    scores = m.calculate()
    assert np.allclose(m.scores, expected, rtol=R_TOLERANCE)


def test_bidir_model():
    S = np.array(
        [
            [
                -1,
                00,
                00,
            ],
            [
                +1,
                -1,
                00,
            ],
            [
                00,
                +1,
                +1,
            ],
            [
                00,
                00,
                -1,
            ],
        ]
    )
    rev = [True for i in range(3)]
    met_names = ["M" + str(x) for x in range(4)]
    expected = np.array([0.175, 0.325, 0.325, 0.175])
    m = Model(S, met_names, rev)
    m.calculate()
    assert np.allclose(m.scores, expected, rtol=R_TOLERANCE)


def test_check_and_reshape_expression_vector():
    S = np.array(
        [
            [-1, 0, 0, 0, 0, 1],
            [1, -1, 0, 0, 0, 0],
            [0, 1, -1, 0, 0, 0],
            [0, 0, 1, -1, 0, 0],
            [0, 0, 0, 1, -1, 0],
            [0, 0, 0, 0, 1, -1],
        ]
    )
    rev = [False for i in range(6)]
    expected = np.ones(6) * 1 / 6
    met_names = ["M" + str(x) for x in range(6)]
    m = Model(S, met_names, rev)
    m.expression_vector = np.ones((2, 3)) * 2
    m._check_and_reshape_expression_vector()
    assert m.expression_vector.shape == (1, 6)


def test_get_subnetworks():
    S = np.array(
        [
            [-1, 0, 0, 0, 0],
            [1, -1, 0, 0, 0],
            [0, 1, -1, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, -1, 1],
            [0, 0, 0, 1, -1],
        ]
    )
    mets = [f"M{i+1}" for i in range(6)]
    rev = [False for i in range(5)]
    model = Model(S, mets, rev)
    subnets = model.get_subnetworks()
    print(test_get_subnetworks)
    assert len(subnets) == 2
    assert len(subnets[0]) == 4
    assert len(subnets[1]) == 2
