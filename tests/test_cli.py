from dataclasses import dataclass
from pathlib import Path
from subprocess import run
from typing import Optional

import pytest
from numpy import allclose, isclose
from pandas import Series, read_csv

from gemcat import cli


@dataclass
class MockNamespace:
    expressionfile: str
    modelfile: str
    expressioncolumn: Optional[str] = None
    baseline: Optional[str] = None
    baselinecolumn: Optional[str] = None
    genefill: Optional[str] = None
    outfile: Optional[str] = None


expression_path = Path("./tests/test_seq/")
model_path = Path("./tests/test_models/")
results_path = Path("./tests/test_results/")


def test_cli_mini():
    args_minimal = MockNamespace(
        expressionfile=str(expression_path / "expression_mini.csv"),
        expressioncolumn=None,
        baseline=None,
        baselinecolumn=None,
        genefill="1.0",
        modelfile=str(model_path / "mini.xml"),
        outfile="./temp_outfile.csv",
    )
    result, outfile = cli.cli_standard(args_minimal)
    assert isinstance(result, Series)
    assert outfile.stem == "temp_outfile"
    assert outfile.suffix == ".csv"
    assert isclose(result["A"], 1.0, atol=0.001)
    assert isclose(result["B"], 0.987, atol=0.001)
    assert len(result) == 4


@pytest.mark.slow
def test_cli_uc_xml():
    args_uc_xml = MockNamespace(
        expressionfile=str(expression_path / "prot_uc_vs_healthy_gid.csv"),
        expressioncolumn="foldchange",
        baseline=None,
        baselinecolumn=None,
        genefill="1.0",
        modelfile=str(model_path / "Recon3D_301.xml"),
        outfile="./temp_outfile.csv",
    )
    result, outfile = cli.cli_standard(args_uc_xml)
    assert isinstance(result, Series)
    assert outfile.stem == "temp_outfile"
    assert outfile.suffix == ".csv"


@pytest.mark.slow
def test_cli_293_xml():
    args_293_xml = MockNamespace(
        expressionfile=str(expression_path / "hek293_rnaseq.csv"),
        expressioncolumn="wtHEK293_1",
        baseline=str(expression_path / "hek293_rnaseq.csv"),
        baselinecolumn="wtHEK293_2",
        genefill="1.0",
        modelfile=str(model_path / "Recon3D_301.xml"),
        outfile="./temp_outfile.csv",
    )
    result, outfile = cli.cli_standard(args_293_xml)
    assert isinstance(result, Series)
    assert len(result) > 1000
    assert outfile.stem == "temp_outfile"
    assert outfile.suffix == ".csv"


@pytest.mark.slow
def test_cli_paper_uc_ns():
    args = MockNamespace(
        expressionfile=str(expression_path / "prot_uc_vs_healthy.csv"),
        expressioncolumn="foldchange",
        baseline=str(expression_path / "prot_uc_vs_healthy.csv"),
        baselinecolumn="base",
        genefill=None,
        modelfile=str(model_path / "Recon3D.json"),
        outfile="./temp_outfile.csv",
    )
    result, outfile = cli.cli_standard(args)
    assert isinstance(result, Series)
    assert len(result) > 1000
    assert outfile.stem == "temp_outfile"
    assert outfile.suffix == ".csv"
    expected = read_csv(results_path / "uc.csv", sep=",", index_col=0).iloc[:, 0]
    assert allclose(expected.values, result.values, atol=0.2)


@pytest.mark.slow
def test_cli_paper_uc_wf():
    args = MockNamespace(
        expressionfile=str(expression_path / "prot_uc_vs_healthy.csv"),
        expressioncolumn="foldchange",
        baseline=str(expression_path / "prot_uc_vs_healthy.csv"),
        baselinecolumn="base",
        genefill=None,
        modelfile=str(model_path / "Recon3D.json"),
        outfile="./temp_outfile.csv",
    )
    result, _ = cli.cli_standard(args)
    assert isinstance(result, Series)
    assert len(result) > 1000
    expected = read_csv(results_path / "uc.csv", sep=",", index_col=0).iloc[:, 0]
    assert allclose(expected.values, result.values, atol=0.2)


@pytest.mark.slow
def test_cli_paper_uc_term():
    out_file = "./temp_outfile.csv"
    run(
        [
            "gemcat",
            str(model_path / "Recon3D.json"),
            str(expression_path / "prot_uc_vs_healthy.csv"),
            "-e",
            "foldchange",
            "-o",
            out_file,
            "-g",
            "1.0",
        ],
        shell=True,
    )
    received = read_csv(out_file, sep=",", index_col=0).iloc[:, 0]
    expected = read_csv(results_path / "uc.csv", sep=",", index_col=0).iloc[:, 0]
    assert allclose(expected.values, received.values, atol=0.2)
    # when run with pytest, maximum deviation is at .09,
    # when run as a normal Python function, it's at 0.01
    # in both cases only the one test function runs
    # weird
