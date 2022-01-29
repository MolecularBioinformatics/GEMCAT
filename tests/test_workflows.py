from pyreporter import workflows as wf
from fixtures import *
import pandas as pd
import numpy as np

regex_tests = 'G\d'

def test_workflow_single(models):
    mini = models['mini']
    genes = [f'G{i}' for i in range(1,6)]
    genes_mapped = pd.Series(dict(zip(genes, [2.]*5)))
    result = wf.workflow_single(mini, genes_mapped)
    assert isinstance(result[0], pd.Series)

def test_workflow_Fang2012(models):
    mini = models['mini']
    genes = [f'G{i}' for i in range(1,6)]
    genes_baseline = pd.Series(dict(zip(genes, [1.]*5)))
    genes_mapped = pd.Series(dict(zip(genes, [2.]*5)))
    result = wf.workflow_Fang2012(
        mini, 
        genes_mapped,
        genes_baseline,
        re_gene = regex_tests
        )
    assert isinstance(result[0], pd.Series)
