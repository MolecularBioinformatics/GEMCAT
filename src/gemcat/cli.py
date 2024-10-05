#!/usr/bin/python

"""
Command line interface functionality
"""

import argparse
import csv
import logging
from copy import deepcopy
from pathlib import Path
from typing import Any, Optional

import cobra
from pandas import DataFrame, Series, read_csv

from .io import load_json_cobra, load_mat_cobra, load_sbml_cobra
from .model_manager import ModelManager
from .workflows import workflow_standard


def not_implemented(whatever: Any):
    """
    Standard function raising a NotImplementedError on demand
    :param whatever: Any input
    :type whatever: Any
    :raises NotImplementedError: Always raises
    """
    raise NotImplementedError()


MODELS = {
    "sbml": load_sbml_cobra,
    "xml": load_sbml_cobra,
    "json": load_json_cobra,
    "csv": not_implemented,
    "mat": load_mat_cobra,
}

model_manager = ModelManager()


def wrong_filetype(any: Any):
    raise NotImplementedError(f"Not implemented for {any}")


def parse_model(model_path: str) -> cobra.Model:
    model_path = Path(model_path)
    throw_for_missing_model(model_path)
    parsing_fn = MODELS[model_path.suffix[1:]]
    return parsing_fn(model_path)


def get_delimiter(file_path, bytes=4096):
    sniffer = csv.Sniffer()
    data = open(file_path, "r").read(bytes)
    delimiter = sniffer.sniff(data).delimiter
    return delimiter


def throw_for_missing_model(model_path: Path):
    """
    Throw FileNotFoundError in case there is no model at the given path.
    :param model_path: _description_
    :type model_path: str
    :raises FileNotFoundError: _description_
    """
    if not model_path.is_file():
        error_str = f"The model file at {model_path} cannot be found."
        logging.error(error_str)
        raise FileNotFoundError(error_str)


def read_expression(expression_file: str) -> DataFrame:
    """
    Parsing the expression file into a DataFrame
    :param expression_file: Path to expression file
    :type expression_file: str
    :raises FileNotFoundError: In case the file does not exist
    :raises ValueError: In tcase the file format is not either CSV or TSV
    :return: DataFrame of the whole expression file
    :rtype: DataFrame
    """
    file_path = Path(expression_file)
    if not file_path.exists():
        error_str = f"File {expression_file} could not be found"
        logging.error(error_str)
        raise FileNotFoundError(error_str)
    if file_path.suffix == ".csv" or file_path.suffix == "tsv":
        delimiter = get_delimiter(file_path)
        content = read_csv(file_path, sep=delimiter, index_col=0)
    else:
        error_str = (
            f"Unsupported file format {file_path.suffix} in file {expression_file}"
        )
        logging.error(error_str)
        raise ValueError(error_str)
    if isinstance(content, Series):
        return DataFrame(content)
    return content


def parse_expression(expression_file: str, col_name: Optional[str]) -> Series:
    """
    Parse file containing expression data into pandas Series
    :param expression_file: Path to expression file
    :type expression_file: str
    :param col_name: ID of the column containing the expression data
    :type col_name: Optional[str]
    :raises ValueError: If expression file cannot be parsed into a series
    :return: _description_
    :rtype: Series
    """
    content = read_expression(expression_file)
    if col_name is None and len(content.columns) > 1:
        raise ValueError(
            """
            If your expression file contains more than 1 column, 
            please provide the name of the column with the expression data
            to the -e flag of the command.
            """
        )
    if col_name is None:
        return content.iloc[:, 0]
    return content.loc[:, col_name]


def get_all_ones(expression: Series) -> Series:
    """
    Return a copy of a pandas series with all entries set to one
    :param expression: Gene expression values
    :type expression: pd.Series
    :return: Identical expression series with all entries set to one
    :rtype: pd.Series
    """
    return Series(index=deepcopy(expression.index), data=1.0)


def parse_outfile(outfile: Optional[str]) -> Path:
    """
    Parse input for output file
    :param outfile: Output file to write to, default to ./results.csv
    :type outfile: Optional[str]
    :return: Path to output file
    :rtype: Path
    """
    if not outfile:
        outfile = "./results.csv"
    outfile = Path(outfile)
    if outfile.suffix not in [".csv", ".tsv"]:
        logging.warning("Cannot handle output format %s .", outfile.suffix)
        outfile = outfile.with_suffix(".csv")
        logging.warning("Saving instead to %s .", outfile)
    return outfile


def build_parser() -> argparse.ArgumentParser:
    """
    Build up CLI parser
    :return: _description_
    :rtype: argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog="gemcat",
        description="GEMCAT tool for metabolomics predictions",
    )

    parser.add_argument(
        "modelfile",
        help=f"Path to model file to use (XML/SBML, JSON, or MAT format), or one of: {model_manager.get_managed_models_str()}",
    )
    parser.add_argument(
        "expressionfile", help="Path to file containing the condition expression data"
    )

    parser.add_argument(
        "-e",
        "--expressioncolumn",
        help="Name of the column containing the condition expression data",
    )
    parser.add_argument(
        "-b", "--baseline", help="File containing expression data for the baseline"
    )
    parser.add_argument(
        "-c",
        "--baselinecolumn",
        help="Name of the column contianing the baseline expression data",
    )
    parser.add_argument(
        "-g",
        "--genefill",
        type=float,
        help="Value to fill in for missing expression values",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Use for more verbose output"
    )
    parser.add_argument("-o", "--outfile", help="Write outputs to the specified file")
    parser.add_argument("-l", "--logfile", help="Write logs to the specified file")

    return parser


def save_to_file(outfile: Path, results: Series) -> None:
    """
    Write results to file
    :param outfile: Path to file to write to
    :type outfile: Path
    :param results: Results to write into the file
    :type results: pd.Series
    :raises ValueError: In case of the file having unsuitable file endings
    """
    if outfile.exists():
        warn_str = f"Overwriting existing file at {outfile}."
        logging.debug(warn_str)
    if outfile.suffix == ".csv":
        results.to_csv(outfile)
    elif outfile.suffix == ".tsv":
        results.to_csv(outfile, sep="\t")
    else:
        # We should never end up here, as we only have .csv and .tsv output
        raise ValueError(
            "An unknown issue occurred with the output file. Output is not being saved."
        )


def cli_standard(args: argparse.Namespace):
    """ "
    CLI interface for the standard workflow
    """
    expression = parse_expression(args.expressionfile, args.expressioncolumn)
    if args.baseline:
        baseline = parse_expression(args.baseline, args.baselinecolumn)
    else:
        print("Empty baseline expression. Defaulting to all ones.")
        baseline = get_all_ones(expression)
    try:
        gene_fill = float(args.genefill)
    except (TypeError, ValueError):
        logging.info("Empty or invalid gene-fill value. Defaulting to 1.0 .")
        gene_fill = 1.0
    if args.modelfile in model_manager.get_managed_models():
        model_path = model_manager.get_model(args.modelfile)
        cobra_model = parse_model(model_path)
    else:
        cobra_model = parse_model(args.modelfile)
    return workflow_standard(
        cobra_model, baseline, expression, gene_fill
    ), parse_outfile(args.outfile)


def main():
    """
    CLI entrypoint and central function
    """
    parser = build_parser()
    args = parser.parse_args()
    verbosity = logging.DEBUG if args.verbose else logging.WARN
    datefmt = "%Y-%m-%d %H:%M:%S"
    fmt = "%(asctime)s %(levelname)s %(message)s"
    if args.logfile:
        logging.basicConfig(
            format=fmt,
            level=verbosity,
            datefmt=datefmt,
            filename=args.logfile,
            encoding="utf-8",
        )
    else:
        logging.basicConfig(
            format=fmt,
            level=verbosity,
            datefmt=datefmt,
        )
    results, outfile = cli_standard(args)
    save_to_file(outfile, results)
