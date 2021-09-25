import pyreporter as pr
import pandas as pd
import numpy as np
import pytest
from fixtures import *

def test_sanity_double_all(models):
    genes_mapped_single = pd.Series({
        'G1': 1,
        'G2': 1,
        'G3': 1,
        'G4': 1,
        'G5': 1,
    })
    genes_mapped_double = pd.Series({
        'G1': 2,
        'G2': 2,
        'G3': 2,
        'G4': 2,
        'G5': 2,
    })
    mini = models['mini']
    result_single = pr.workflows.workflow_single(mini, genes_mapped_single)
    result_double = pr.workflows.workflow_single(mini, genes_mapped_double)
    assert np.allclose(result_single[0].values, result_double[0].values)