from pyreporter import workflows as wf
from fixtures import *
import pandas as pd
import numpy as np

def test_workflow_single(models):
    mini = models['mini']
    genes_mapped = pd.Series({
        'G1': 2,
        'G2': 2,
        'G3': 2,
        'G4': 2,
        'G5': 2,
    })
    result = wf.workflow_single(mini, genes_mapped)
    assert isinstance(result[0], pd.Series)

