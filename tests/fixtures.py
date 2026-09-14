from pathlib import Path

import cobra
import numpy as np
import pytest
import scipy.sparse as sp


def to_dense(matrix):
    """
    Densify a matrix for comparison against a dense expectation.

    Use this at every assertion site rather than comparing directly. GEMCAT returns sparse
    matrices, and `np.allclose(sparse, dense)` either raises NotImplementedError or -- worse,
    in some shapes -- silently densifies and passes for the wrong reason.
    """
    return matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)


@pytest.fixture
def models():
    modelpath = Path("./tests/test_models/")
    modelpaths = {
        "mini": modelpath / "mini.xml",
        "mini_redox": modelpath / "mini_redox.xml",
        "mini_reversible": modelpath / "mini_reversible.xml",
        "mini_complex_gpr": modelpath / "mini_complex_gpr.xml",
        "mini_complex_gpr_gids": modelpath / "mini_complex_gpr_gid_version.xml",
    }
    models = {
        name: cobra.io.read_sbml_model(p.as_posix()) for name, p in modelpaths.items()
    }
    return models


@pytest.fixture
def model_files_sbml():
    modelpath = Path("./tests/test_models/")
    modelpaths = {
        "mini": modelpath / "mini.xml",
        "mini_redox": modelpath / "mini_redox.xml",
        "mini_reversible": modelpath / "mini_reversible.xml",
        "mini_complex_gpr": modelpath / "mini_complex_gpr.xml",
        "mini_complex_gpr_gids": modelpath / "mini_complex_gpr_gid_version.xml",
    }
    modelpaths = {name: Path(loc) for name, loc in modelpaths.items()}
    return modelpaths


@pytest.fixture
def model_files_json():
    modelpath = Path("./tests/test_models/")
    modelpaths = {
        "mini": modelpath / "mini.json",
        "mini_redox": modelpath / "mini_redox.json",
        "mini_reversible": modelpath / "mini_reversible.json",
        "mini_complex_gpr": modelpath / "mini_complex_gpr.json",
        "mini_complex_gpr_gids": modelpath / "mini_complex_gpr_gid_version.json",
    }
    modelpaths = {name: Path(loc) for name, loc in modelpaths.items()}
    return modelpaths


@pytest.fixture
def S_examples():
    """
    Stoichiometric matrices used as transform *input*, so they are returned sparse.

    The literals stay dense and hand-written because they document the maths; only the
    returned values are converted. Use S_examples_dense if you need the raw array.
    """
    return {
        name: sp.csr_array(value.astype(float))
        for name, value in _S_examples_dense().items()
    }


@pytest.fixture
def S_examples_dense():
    return _S_examples_dense()


def _S_examples_dense():
    S = {}
    S["complex"] = np.array(
        [
            [-1, +1, -1, +1, -1, -1, +1, +1, -1, 00, 00, 00, 00, 00, 00, 00, 00, 00],
            [+1, -1, 00, 00, 00, 00, 00, 00, 00, +1, +1, 00, 00, 00, 00, 00, 00, 00],
            [00, 00, +1, -1, 00, 00, 00, 00, 00, -1, 00, +1, +1, 00, 00, 00, 00, 00],
            [00, 00, 00, 00, +1, 00, 00, 00, 00, 00, -1, -1, 00, -1, +1, 00, 00, 00],
            [00, 00, 00, 00, 00, +1, -1, 00, 00, 00, 00, 00, -1, +1, -1, -1, +1, +1],
            [00, 00, 00, 00, 00, 00, 00, -1, 00, 00, 00, 00, 00, 00, 00, +1, -1, 00],
            [00, 00, 00, 00, 00, 00, 00, 00, +1, 00, 00, 00, 00, 00, 00, 00, 00, -1],
        ]
    )
    S["circular"] = np.array(
        [
            [-1, 0, 0, 0, 0, 1],
            [1, -1, 0, 0, 0, 0],
            [0, 1, -1, 0, 0, 0],
            [0, 0, 1, -1, 0, 0],
            [0, 0, 0, 1, -1, 0],
            [0, 0, 0, 0, 1, -1],
        ]
    )
    S["linear"] = np.array(
        [
            [-1, 0, 0, 0, 0, 0],
            [1, -1, 0, 0, 0, 0],
            [0, 1, -1, 0, 0, 0],
            [0, 0, 1, -1, 0, 0],
            [0, 0, 0, 1, -1, 0],
            [0, 0, 0, 0, 1, 0],
        ]
    )
    S["bidir_linear"] = np.array(
        [
            [
                -1,
                +1,
                00,
                00,
                00,
                00,
            ],
            [
                +1,
                -1,
                -1,
                +1,
                00,
                00,
            ],
            [
                00,
                00,
                +1,
                -1,
                -1,
                +1,
            ],
            [
                00,
                00,
                00,
                00,
                +1,
                -1,
            ],
        ]
    )
    return S


@pytest.fixture
def A_examples():
    A = {}
    A["complex"] = np.array(
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
    A["circular"] = np.array(
        [
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0, 0],
        ]
    )
    A["linear"] = np.array(
        [
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 0],
        ]
    )
    A["bidir_linear"] = np.array(
        [
            [000, 1 / 1, 000, 000],
            [1 / 2, 000, 1 / 2, 000],
            [000, 1 / 2, 000, 1 / 2],
            [000, 000, 1 / 1, 000],
        ]
    )
    return A


@pytest.fixture
def S_models():
    S_models = {}
    S_models["mini"] = np.array(
        [
            [-1, -1, 00, 00],
            [+1, 00, -1, 00],
            [00, +1, 00, -1],
            [00, 00, +1, +1],
        ]
    )
    return S_models
