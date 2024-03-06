from dataclasses import dataclass
from typing import Optional

from numpy import isclose
from pandas import Series

from gemcat import cli
from pathlib import Path


@dataclass
class MockNamespace:
    expressionfile: str
    modelfile: str
    expressioncolumn: Optional[str] = None
    baseline: Optional[str] = None
    baselinecolumn: Optional[str] = None
    genefill: Optional[str] = None
    ranking: Optional[str] = None
    adjacency: Optional[str] = None
    outfile: Optional[str] = None
    integration: Optional[str] = None

expression_path = Path('./tests/test_seq/')
model_path = Path('./tests/test_models/')

args_minimal = MockNamespace(
    expressionfile=str(expression_path / "expression_mini.csv"),
    expressioncolumn=None,
    baseline=None,
    baselinecolumn=None,
    genefill="1.0",
    modelfile=str(model_path / "mini.xml"),
    ranking=None,
    adjacency=None,
    outfile="./temp_outfile.csv",
    integration="means",
)

args_uc_xml = MockNamespace(
    expressionfile=str(expression_path / "test_prot_uc_vs_healthy.csv"),
    expressioncolumn="foldchange",
    baseline=None,
    baselinecolumn=None,
    genefill="1.0",
    modelfile=str(model_path / "recon3d.xml"),
    ranking=None,
    adjacency=None,
    outfile="./temp_outfile.csv",
    integration="means",
)


def test_cli_mini():
    result, outfile = cli.cli_standard(args_minimal)
    print(result)
    assert isinstance(result, Series)
    assert outfile.stem == "temp_outfile"
    assert outfile.suffix == ".csv"
    assert isclose(result["A"], 1.0, atol=0.001)
    assert isclose(result["B"], 0.987, atol=0.001)
    assert len(result) == 4


def test_cli_uc_xml():
    result, outfile = cli.cli_standard(args_uc_xml)
    print(result)
    assert isinstance(result, Series)
    assert outfile.stem == "temp_outfile"
    assert outfile.suffix == ".csv"
