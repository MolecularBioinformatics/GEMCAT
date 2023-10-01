import argparse
import logging
from pathlib import Path
from typing import Optional

import cobra
from pandas import DataFrame, Series, read_csv

from .AdjacencyTransformation import AdjacencyTransformation, ATPureAdjacency
from .io import load_sbml_cobra
from .Model import Model
from .PageRank import PageRankNX, Ranking
from .workflows import workflow_Fang2012

ADJACENCIES = {
    "pure": ATPureAdjacency,
}
ALLOWED_ADJ = ", ".join(ADJACENCIES.keys())

RANKINGS = {
    "pagerank": PageRankNX,
}
ALLOWED_RANKINGS = ", ".join(RANKINGS.keys())


def not_implemented(whatever):
    raise NotImplementedError()


MODELS = {
    ".sbml": load_sbml_cobra,
    ".xml": load_sbml_cobra,
    ".json": not_implemented,
    ".csv": not_implemented,
    ".mat": not_implemented,
}


def parse_cobra_model(model_path: str) -> cobra.Model:
    model_path = Path(model_path)
    if not model_path.exists():
        error_str = f"The model file at {model_path} cannot be found."
        logging.error(error_str)
        raise FileNotFoundError(error_str)
    return cobra.io.read_sbml_model(model_path.as_posix())


def read_expression(expression_file: str) -> DataFrame:
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
    content = read_expression(expression_file)
    if not content.shape[1] or (content.shape[1] == 1):
        content = content.squeeze()
    else:
        content = content.loc[:, col_name]
    if not isinstance(content, Series):
        error_str = "Expression file must be in the format of a pandas Series, mapping a gene column to a value column."
        logging.error(error_str)
        raise ValueError(error_str)
    return content


def parse_adjacency(adjacency: Optional[str]) -> AdjacencyTransformation:
    if adjacency is None:
        return ATPureAdjacency
    try:
        return ADJACENCIES[adjacency]
    except KeyError:
        error_str = (
            f"Adjacency model {adjacency} does not exist. Allowed models: {ALLOWED_ADJ}"
        )
        logging.error(error_str)
        raise ValueError(error_str)


def parse_ranking(ranking: Optional[str]) -> Ranking:
    if ranking is None:
        return PageRankNX
    try:
        return RANKINGS[ranking]
    except KeyError:
        error_str = f"{ranking} is not a valid ranking method. Available methods are: {ALLOWED_RANKINGS}"
        logging.error(error_str)
        raise ValueError(error_str)


def get_all_ones(expression: Series) -> Series:
    return Series(index=expression.index, data=1.0)


def parse_outfile(outfile: Optional[str]) -> Path:
    if not outfile:
        outfile = "./results.csv"
    outfile = Path(outfile)
    if not outfile.suffix in [".csv", ".tsv"]:
        logging.warning(f"Cannot handle output format {outfile.suffix} .")
        outfile = outfile.with_suffix(".csv")
        logging.warning(f"Saving instead to {outfile} .")
    return outfile


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="PyReporter",
        description="PyReporter tool for metabolomics predictions",
    )
    parser.add_argument("expressionfile")
    parser.add_argument("modelfile")

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


def cli_fang2012(args: argparse.Namespace):
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

    model = parse_cobra_model(args.modelfile)
    adjacency = parse_adjacency(args.adjacency)
    ranking = parse_ranking(args.ranking)
    outfile = parse_outfile(args.outfile)

    results = workflow_Fang2012(
        model,
        baseline,
        expression,
        adjacency=adjacency,
        ranking=ranking,
        gene_fill=gene_fill,
    )
    return results, outfile


def main():
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
    results, outfile = cli_fang2012(args)
    save_to_file(outfile, results)
