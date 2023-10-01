import numpy as np
import pandas as pd
from fixtures import *

from prankme import workflows as wf


def test_workflow_single(models):
    mini = models["mini"]
    genes = [f"G{i}" for i in range(1, 6)]
    genes_mapped = pd.Series(dict(zip(genes, [2.0] * 5)))
    result = wf.workflow_single(mini, genes_mapped)
    assert isinstance(result, pd.Series)


def test_workflow_Fang2012(models):
    mini = models["mini"]
    genes = [f"G{i}" for i in range(1, 6)]
    genes_baseline = pd.Series(dict(zip(genes, [1.0] * 5)))
    genes_mapped = pd.Series(dict(zip(genes, [2.0] * 5)))
    result = wf.workflow_Fang2012(
        mini,
        genes_mapped,
        genes_baseline,
    )
    assert isinstance(result, pd.Series)
