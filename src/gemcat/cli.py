#!/usr/bin/python

"""
Command line interface functionality
"""

import argparse
import logging
from pathlib import Path
from typing import Any, Optional

import cobra
from pandas import DataFrame, Series, read_csv

from .adjacency_transformation import AdjacencyTransformation, ATPureAdjacency
from .expression import (
    ExpressionIntegration,
    ExpressionMapSingleAverage,
    GeometricAndAverageMeans,
    read_gpr_strings_from_cobra,
)
from .io import convert_cobra_model, load_sbml_cobra
from .ranking import PagerankNX, Ranking

ADJACENCIES = {
    "pure": ATPureAdjacency,
}
ALLOWED_ADJ = ", ".join(ADJACENCIES.keys())

RANKINGS = {
    "Pagerank": PagerankNX,
}

EXPRESSION_INTEGRATIONS = {
    "means": GeometricAndAverageMeans,
    "average": ExpressionMapSingleAverage,
}

ALLOWED_RANKINGS = ", ".join(RANKINGS.keys())


def not_implemented(whatever: Any):
    """
    Standard function raising a NotImplementedError on demand
    :param whatever: Any input
    :type whatever: Any
    :raises NotImplementedError: Always raises
    """
    raise NotImplementedError()


MODELS = {
    ".sbml": load_sbml_cobra,
    ".xml": load_sbml_cobra,
    ".json": not_implemented,
    ".csv": not_implemented,
    ".mat": not_implemented,
}


def parse_cobra_model(model_path: str) -> cobra.Model:
    """
    Parse selected cobra model
    :param model_path: Path to model file
    :type model_path: str
    :raises FileNotFoundError: If model file does not exist at path
    :return: Loaded cobra model
    :rtype: cobra.Model
    """
    model_path = Path(model_path)
    if not model_path.exists():
        error_str = f"The model file at {model_path} cannot be found."
        logging.error(error_str)
        raise FileNotFoundError(error_str)
    return cobra.io.read_sbml_model(model_path.as_posix())


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
    if file_path.suffix == ".csv":
        content = read_csv(file_path, sep=",", index_col=0)
    elif file_path.suffix == ".tsv":
        content = read_csv(file_path, sep="\t", index_col=0)
    else:
        error_str = f"Unknown file format {file_path.suffix} in file {expression_file}"
        logging.error(error_str)
        raise ValueError(error_str)
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
    if not content.shape[1] or (content.shape[1] == 1):
        content = content.squeeze()
    else:
        content = content.loc[:, col_name]
    if not isinstance(content, Series):
        error_str = """
        Expression file must be in the format of a pandas Series, 
        mapping a gene column to a value column.
        """
        logging.error(error_str)
        raise ValueError(error_str)
    return content


def parse_integration(integration: str) -> ExpressionIntegration:
    """
    Choice function for expression integration algorithms
    :param integration: Name for expression integration algorithm
    :type integration: str
    :raises ValueError: In case of invalid name
    :return: Selected expression integration algorithm
    :rtype: ExpressionIntegration
    """
    if integration is None:
        return GeometricAndAverageMeans
    try:
        return EXPRESSION_INTEGRATIONS[integration]
    except KeyError as original_err:
        error_str = f"""
        Expression integration method {integration} does not exist. 
        Allowed methods: {EXPRESSION_INTEGRATIONS}
        """
        logging.error(error_str)
        raise ValueError(error_str) from original_err


def parse_adjacency(adjacency: Optional[str]) -> AdjacencyTransformation:
    """
    Choice function for adjacency calculation algorithms
    :param adjacency: Name for adjacency calculation algorithm
    :type adjacency: Optional[str]
    :raises ValueError: In case of invalid name
    :return: Selected adjacency transformation algorithm
    :rtype: AdjacencyTransformation
    """
    if adjacency is None:
        return ATPureAdjacency
    try:
        return ADJACENCIES[adjacency]
    except KeyError as original_err:
        error_str = f"""
        Adjacency model {adjacency} does not exist. 
        Allowed models: {ALLOWED_ADJ}
        """
        logging.error(error_str)
        raise ValueError(error_str) from original_err


def parse_ranking(ranking: Optional[str]) -> Ranking:
    """
    Choice function for ranking algorithms
    :param ranking: Name of the ranking algorithm to use
    :type ranking: Optional[str]
    :raises ValueError: In case of unknown ranking algorithm name
    :return: Ranking algorithm
    :rtype: Ranking
    """
    if ranking is None:
        return PagerankNX
    try:
        return RANKINGS[ranking]
    except KeyError as original_err:
        error_str = f"""
        {ranking} is not a valid ranking method. 
        Available methods are: {ALLOWED_RANKINGS}
        """
        logging.error(error_str)
        raise ValueError(error_str) from original_err


def get_all_ones(expression: Series) -> Series:
    """
    Return a copy of a pandas series with all entries set to one
    :param expression: Gene expression values
    :type expression: pd.Series
    :return: Identical expression series with all entries set to one
    :rtype: pd.Series
    """
    return Series(index=expression.index, data=1.0)


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
    parser.add_argument("expressionfile")
    parser.add_argument("modelfile")

    parser.add_argument("-i", "--integration")
    parser.add_argument("-e", "--expressioncolumn")
    parser.add_argument("-r", "--ranking")
    parser.add_argument("-b", "--baseline")
    parser.add_argument("-c", "--baselinecolumn")
    parser.add_argument("-a", "--adjacency")
    parser.add_argument("-g", "--genefill")
    parser.add_argument("-v", "--verbose")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-l", "--logfile")

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
            "An unknown issue occured with the output file. Output is not being saved."
        )


def cli_standard(args: argparse.Namespace):
    """ "
    CLI interface for the standard workflow
    """
    expression = parse_expression(args.expressionfile, args.expressioncolumn)
    if args.baseline:
        baseline = parse_expression(args.baseline, args.baseline_column)
    else:
        print("Empty baseline expression. Defaulting to all ones.")
        baseline = get_all_ones(expression)
    try:
        gene_fill = float(args.genefill)
    except (TypeError, ValueError):
        print("Empty or invalid gene-fill value. Defaulting to 1.0 .")
        gene_fill = 1.0

    cobra_model = parse_cobra_model(args.modelfile)
    model = convert_cobra_model(cobra_model)
    model.adjacency_transformation = parse_adjacency(args.adjacency)
    model.ranking = parse_ranking(args.ranking)
    gpr, rxn_gene_mapping = read_gpr_strings_from_cobra(cobra_model)
    integration = parse_integration(args.integration)
    expression_obj = integration(gpr, rxn_gene_mapping, expression, gene_fill=gene_fill)
    baseline_obj = integration(gpr, rxn_gene_mapping, baseline, gene_fill=gene_fill)
    model.load_expression(expression_obj)
    results = model.calculate()
    model.load_expression(baseline_obj)
    results_baseline = model.calculate()

    return results / results_baseline, parse_outfile(args.outfile)


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
