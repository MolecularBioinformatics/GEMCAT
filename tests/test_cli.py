from prankme import cli
from typing import Optional
from dataclasses import dataclass
from pandas import Series
from numpy import isclose


@dataclass
class MockNamespace:
    expressionfile: Optional[str] = None
    expressioncolumn: Optional[str] = None
    baseline: Optional[str] = None
    baselinecolumn: Optional[str] = None
    genefill: Optional[str] = None
    modelfile: Optional[str] = None
    ranking: Optional[str] = None
    adjacency: Optional[str] = None
    outfile: Optional[str] = None


args_minimal = MockNamespace(
    expressionfile="./tests/test_seq/expression_mini.csv",
    expressioncolumn=None,
    baseline=None,
    baselinecolumn=None,
    genefill="1.0",
    modelfile="./tests/test_models/mini.xml",
    ranking=None,
    adjacency=None,
    outfile="./temp_outfile.csv",
)


def test_cli_workflow():
    result, outfile = cli.cli_fang2012(args_minimal)
    print(result)
    assert isinstance(result, Series)
    assert outfile.stem == "temp_outfile"
    assert outfile.suffix == ".csv"
    assert isclose(result["A"], 1.0, atol=0.001)
    assert isclose(result["B"], 0.987, atol=0.001)
    assert len(result) == 4
