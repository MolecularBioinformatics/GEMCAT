import numpy as np
import sparse

import gemcat.ranking as pr
from gemcat.utils import all_close

TOLERANCE = 10**-3


def test_simple_linear():
    A = sparse.COO(
        [
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 0],
        ]
    )
    result = pr.PagerankNX.simple_pagerank_nx(A)
    expected = [0.061, 0.112, 0.156, 0.193, 0.225, 0.252]
    assert all_close(result, expected, atol=TOLERANCE)


def test_simple_circular_nx():
    A = sparse.COO(
        [
            [
                0,
                1,
                0,
                0,
                0,
            ],
            [
                0,
                0,
                1,
                0,
                0,
            ],
            [
                0,
                0,
                0,
                1,
                0,
            ],
            [
                0,
                0,
                0,
                0,
                1,
            ],
            [
                1,
                0,
                0,
                0,
                0,
            ],
        ]
    )
    result = pr.PagerankNX.simple_pagerank_nx(A)
    expected = np.ones(5) * 0.2
    assert all_close(result, expected, atol=TOLERANCE)


def test_simple_bidir():
    A = sparse.COO(
        [
            [0.0, 1.0, 0.0, 0.0],
            [0.5, 0.0, 0.5, 0.0],
            [0.0, 0.5, 0.0, 0.5],
            [0.0, 0.0, 1.0, 0.0],
        ]
    )
    expected = sparse.COO([0.175, 0.325, 0.325, 0.175])
    result = pr.PagerankNX.simple_pagerank_nx(A)
    assert all_close(result, expected, atol=TOLERANCE)


def test_complex():
    A = sparse.COO(
        [
            [
                000,
                1 / 5,
                1 / 5,
                1 / 5,
                1 / 5,
                000,
                1 / 5,
            ],
            [
                1 / 1,
                000,
                000,
                000,
                000,
                000,
                000,
            ],
            [
                1 / 2,
                1 / 2,
                000,
                000,
                000,
                000,
                000,
            ],
            [
                000,
                1 / 3,
                1 / 3,
                000,
                1 / 3,
                000,
                000,
            ],
            [
                1 / 4,
                000,
                1 / 4,
                1 / 4,
                000,
                1 / 4,
                000,
            ],
            [
                1 / 2,
                000,
                000,
                000,
                1 / 2,
                000,
                000,
            ],
            [
                000,
                000,
                000,
                000,
                1 / 1,
                000,
                000,
            ],
        ]
    )
    expected = [0.280, 0.159, 0.139, 0.108, 0.184, 0.061, 0.069]
    result = pr.PagerankNX.simple_pagerank_nx(A)
    assert all_close(result, expected, atol=TOLERANCE)
