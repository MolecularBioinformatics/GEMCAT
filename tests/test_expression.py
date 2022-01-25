import pyreporter.Expression as ex
import pandas as pd
import numpy as np
from fixtures import *

R_TOLERANCE = 10**-5
gpr_mini = {
    'R1': ['G1'],
    'R2': ['G2'],
    'R3': ['G3'],
    'R4': ['G4'],
}
genes = [f'G{i}' for i in range(1, 9)]
gene_vals = list(range(1,9))
genes_mini_complex = pd.Series(dict(zip(genes, gene_vals)))

def test_read_simple_gpr_from_cobra(models):
    gpr = ex.read_simple_gpr_from_cobra(models['mini'])
    assert gpr == gpr_mini

def test_ExpressionMapSingleAverage(models):
    gene_map = pd.Series({'G1': 2, 'G2': 2, 'G3': 2, 'G4': 2})
    gpr = ex.read_simple_gpr_from_cobra(models['mini'])
    exp = ex.ExpressionMapSingleAverage(gene_map, gpr)
    expected = np.array([2, 2, 2, 2])
    assert np.allclose(exp.get_mapped_values(), expected)

def test_ExpressionMapFang2012_simple(models):
    exp = ex.ExpressionFang2012(
        models['mini_complex_gpr'], 
        genes_mini_complex,
        'G\d'
    )
    assert isinstance(exp.gpr, pd.Series) 
    assert np.isclose(exp.mapped_values['R4'], 4., rtol=R_TOLERANCE)

def test_ExpressionMapFang2012_complex(models):
    exp = ex.ExpressionFang2012(
        models['mini_complex_gpr'], 
        genes_mini_complex,
        'G\d'
    )
    assert isinstance(exp.gpr, pd.Series) 
    assert np.isclose(exp.mapped_values['R3'], 15.4157, rtol=R_TOLERANCE)

def test_ExpressionModifiedFang2012Single_complex(models):
    exp = ex.ExpressionModifiedFang2012Single(
        models['mini_complex_gpr'], 
        genes_mini_complex,
        'G\d'
    )
    assert isinstance(exp.gpr, pd.Series) 
    assert np.isclose(exp.mapped_values['R3'], 15.41574, rtol=R_TOLERANCE)