from pathlib import Path

import numpy as np
import pytest
from pandas import Series

import gemcat as pr

base_path = Path("./tests")
model_path = base_path / "test_models"
seeds_model = model_path / "seeds_test.csv"
mini_path = model_path / "mini.xml"
expression_path = "./tests/test_seq/expression_mini.csv"
model_seeds = model_path / "seeds_test.csv"

eps = 10**-5


def test_seed_load():
    seeds = [1.0, 2.0, 3.0, 4.0]
    model = pr.io.load_sbml_cobra(mini_path)
    model.load_metabolite_seeds(seeds)
    assert model.seeds == seeds


def test_wrong_seed_load():
    seeds = [1.0, 2.0, 3.0, 4.0, 5.0]
    model = pr.io.load_sbml_cobra(mini_path)
    with pytest.raises(ValueError):
        model.load_metabolite_seeds(seeds)


# found here: https://www.briggsby.com/personalized-Pagerank
def test_seed_integration():
    expected = Series(
        {
            "A": 0.226558,
            "B": 0.218752,
            "C": 0.432226,
            "D": 0.122464,
        }
    )
    model = pr.io.load_csv(seeds_model)
    seeds = [0.0, 0.0, 1.0, 0.0]
    gene_lvls = Series(dict(zip(["G1", "G2", "G3", "G4", "G5"], 5 * [1.0])))
    gpr = {
        "R1": ["G1"],
        "R2": ["G2"],
        "R3": ["G2"],
        "R4": ["G3"],
        "R5": ["G3"],
        "R6": ["G4"],
        "R7": ["G5"],
    }
    ex = pr.expression.ExpressionMapSingleAverage(gene_lvls, gpr)
    model.load_expression(ex)
    model.load_metabolite_seeds(seeds)
    results = model.calculate()
    assert isinstance(results, Series)
    assert np.allclose(results.values, expected.values)


def test_no_seed():
    expected = Series(
        {
            "A": 0.257,
            "B": 0.248,
            "C": 0.357,
            "D": 0.139,
        }
    )
    model = pr.io.load_csv(seeds_model)
    gene_lvls = Series(dict(zip(["G1", "G2", "G3", "G4", "G5"], 5 * [1.0])))
    gpr = {
        "R1": ["G1"],
        "R2": ["G2"],
        "R3": ["G2"],
        "R4": ["G3"],
        "R5": ["G3"],
        "R6": ["G4"],
        "R7": ["G5"],
    }
    ex = pr.expression.ExpressionMapSingleAverage(gene_lvls, gpr)
    model.load_expression(ex)
    results = model.calculate()
    assert model.seeds is None
    assert np.allclose(results.values, expected.values, atol=10**-2)
