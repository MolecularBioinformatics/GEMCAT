import pytest
import pyreporter as pr

model_path = './tests/test_models/mini.xml'
expression_path = './tests/test_seq/expression_mini.csv'


def test_seed_load():
    seeds = [1., 2., 3., 4.]
    model = pr.io.load_sbml_cobra(model_path)
    model.load_metabolite_seeds(seeds)
    assert (model.seeds == seeds).all()

def test_wrong_seed_load():
    seeds = [1., 2., 3., 4., 5.]
    model = pr.io.load_sbml_cobra(model_path)
    with pytest.raises(ValueError):
        model.load_metabolite_seeds(seeds)