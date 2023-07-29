import argparse
from cobra import Model as cobraModel
from pathlib import Path
from pandas import Series, read_csv
from typing import Optional
from warnings import warn

from pyreporter.workflows import workflow_Fang2012
from pyreporter.AdjacencyTransformation import AdjacencyTransformation, ATPureAdjacency
from pyreporter.io import load_sbml_cobra, convert_cobra_model
from pyreporter.Model import Model
from pyreporter.PageRank import PageRankNX, Ranking


ADJACENCIES = {
    'pure': ATPureAdjacency,
}
ALLOWED_ADJ = ', '.join(ADJACENCIES.keys())

RANKINGS = {
    'pagerank': PageRankNX,
}
ALLOWED_RANKINGS = ', '.join(RANKINGS.keys())

def not_implemented(whatever):
    raise NotImplementedError()

MODELS = {
    '.sbml': load_sbml_cobra,
    '.xml': load_sbml_cobra, 
    '.json': not_implemented,
    '.csv': not_implemented,
    '.mat': not_implemented,
}

def parse_model(model_path: str) -> cobraModel:
    model_file = Path(model_path)
    if not model_file.exists():
        raise FileNotFoundError('The model file cannot be found.')
    try:
        fn = MODELS[model_file.suffix]
    except KeyError:
        raise ValueError(f'Invalid file type {model_file.suffix}')
    return fn(model_file)

def parse_expression(expression_file: str) -> Series:
    file_path = Path(expression_file)
    if not file_path.exists():
        raise FileNotFoundError(f'File {expression_file} could not be found')
    if file_path.suffix == '.csv':
        content = read_csv(file_path, sep=',', index_col=0)
    elif file_path.suffix == '.tsv':
        content = read_csv(file_path, sep='\t', index_col=0)
    else:
        raise ValueError(f'Unknown file format {file_path.suffix} in file {expression_file}')
    content = content.squeeze()
    if not isinstance(content, Series):
        raise ValueError('Expression file must be in the format of a pandas Series, mapping a gene column to a value column.')
    return content

def parse_adjacency(adjacency: Optional[str]) -> AdjacencyTransformation:
    if adjacency is None:
        return ATPureAdjacency
    try:
        return ADJACENCIES[adjacency]
    except KeyError:
        raise ValueError(f'Adjacency model {adjacency} does not exist. Allowed models: {ALLOWED_ADJ}')

def parse_ranking(ranking: Optional[str]) -> Ranking:
    if ranking is None:
        return PageRankNX
    try:
        return RANKINGS[ranking]
    except KeyError:
        raise ValueError(f'{ranking} is not a valid ranking method. Available methods are: {ALLOWED_RANKINGS}')

def get_all_ones(expression: Series) -> Series:
    return Series(index=expression.index, data = 1.)

def parse_outfile(outfile: Optional[str]) -> Path:
    outfile = Path(outfile)
    if not outfile:
        outfile = Path('./results.csv')
    if not outfile.suffix in ['.csv', '.tsv']:
        warn(f'Cannot handle output format {outfile.suffix} .')
        outfile = outfile.with_suffix('.csv')
        warn(f'Saving instead to {outfile} .')
    return outfile

def main():
    parser = argparse.ArgumentParser(
        prog='PyReporter',
        description='PyReporter tool for metabolomics predictions',
    )
    parser.add_argument('expression_file')
    parser.add_argument('model_file')

    parser.add_argument('-r', '--ranking')
    parser.add_argument('-b', '--baseline')
    parser.add_argument('-a', '--adjacency')
    parser.add_argument('-g', '--genefill')
    parser.add_argument('-v', '--verbose')
    parser.add_argument('-o', '--outfile')

    args = parser.parse_args()
    
    expression = parse_expression(args.expression_file)
    if args.baseline:
       baseline = parse_expression(args.baseline)
    else:
        print('Empty baseline expression. Defaulting to all ones.')
        baseline = get_all_ones(expression)
    
    try:
        gene_fill = float(args.genefill)
    except (TypeError, ValueError):
        print('Empty or invalid gene-fill value. Defaulting to 1.0 .')
        gene_fill = 1.

    model = parse_model(args.model_file)
    adjacency = parse_adjacency(args.adjacency)
    ranking = parse_ranking(args.ranking)
    outfile = parse_outfile(args.outfile)
    
    results = workflow_Fang2012(
        model,
        baseline,
        expression,
        adjacency = adjacency,
        ranking = ranking,
        gene_fill = gene_fill,
    )
    
    if outfile.suffix == 'csv':
        results.to_csv(outfile)
    elif outfile.suffix == 'tsv':
        results.to_csv(outfile, sep='\t')
