import argparse
from pathlib import Path
from pandas import Series, read_csv

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
    'sbml': load_sbml_cobra, 
    'json': not_implemented, 
    'csv': not_implemented,
}

def parse_model(model_path: str) -> Model:
    model_file = Path(model_path)
    if not model_file.exists():
        raise FileNotFoundError('The model file cannot be found.')
    try:
        fn = MODELS[model_file.suffix]
    except KeyError:
        raise ValueError(f'Invalid file type .{model_file.suffix}')
    native_model = fn(model_file)
    return convert_cobra_model(native_model)

def parse_adjacency(adjacency: str) -> AdjacencyTransformation:
    try:
        return ADJACENCIES[adjacency]
    except KeyError:
        raise ValueError(f'Adjacency model {adjacency} does not exist. Allowed models: {ALLOWED_ADJ}')

def parse_expression(expression_file: str) -> Series:
    file_path = Path(expression_file)
    if not file_path.exists():
        raise FileNotFoundError(f'File {expression_file} could not be found')
    if file_path.suffix == 'csv':
        content = read_csv(file_path, sep=',')
    elif file_path.suffix == 'tsv':
        content = read_csv(file_path, sep='\t')
    else:
        raise ValueError('Unknown file format {file_path.suffix} in file {expression_file}')
    if not isinstance(content, Series):
        raise ValueError('Expression file must be in the format of a pandas Series, mapping a gene column to a value column.')
    return content

def parse_ranking(ranking: str) -> Ranking:
    try:
        return RANKINGS[ranking]
    except KeyError:
        raise ValueError(f'{ranking} is not a valid ranking method. Available methods are: {ALLOWED_RANKINGS}')

def get_all_ones(expression: Series):
    return Series(index=expression.index, data = 1.)

def main():
    parser = argparse.ArgumentParser(
        prog='PyReporter',
        description='PyReporter tool for metabolomics predictions',
    )
    parser.add_argument('expression_file')
    parser.add_argument('model_file')

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

    outfile = args.outfile
    if not outfile:
        outfile = Path('./results.csv')
    if not outfile.suffix in ['csv', 'tsv']:
        print(f'Cannot handle output format {outfile.suffix} .')
        outfile = outfile.with_suffix('.csv')
        print('Saving instead to {outfile} .')
    
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
