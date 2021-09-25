import pytest
import pyreporter.PageRank as pr
import numpy as np

TOLERANCE = 10 ** -3

def test_simple_linear():
    A = np.array([
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0],
    ])
    result = pr.run_PR_nx(A)
    expected = [.061, .112, .156, .193, .225, .252]
    assert np.allclose(result, expected, atol=TOLERANCE)

def test_simple_circular_nx():
    A = np.array([
        [0, 0, 1, 0, 0, ],
        [1, 0, 0, 0, 0, ],
        [0, 1, 0, 0, 0, ],
        [0, 0, 0, 0, 1, ],
        [0, 0, 0, 1, 0, ],
    ])
    result = pr.run_PR_nx(A)
    expected = np.ones(5) * .2
    assert np.allclose(result, expected, atol=TOLERANCE)

def test_simple_bidir():
    A = np.array([
        [0. , 1. , 0. , 0. ],
        [0.5, 0. , 0.5, 0. ],
        [0. , 0.5, 0. , 0.5],
        [0. , 0. , 1. , 0. ],
    ])
    expected = np.array([.175, .325, .325, .175])
    result = pr.run_PR_nx(A)
    assert np.allclose(result, expected, atol=TOLERANCE) 

# def test_prs_algos():
#     algos = [
#         'piteration', 'diteration', 'lanczos',
#         'bicgstab', 'RH',
#     ]
#     A = np.array([
#         [0, 0, 1, 0, 0, ],
#         [1, 0, 0, 0, 0, ],
#         [0, 1, 0, 0, 0, ],
#         [0, 0, 0, 0, 1, ],
#         [0, 0, 0, 1, 0, ],
#     ])
#     expected = np.ones(5) * .2
#     for algo in algos:
#         result = prs(A, solver=algo)
#         assert np.allclose(result, expected, atol=TOLERANCE)
