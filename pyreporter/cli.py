from pyreporter.workflows import workflow_Fang2012
from pyreporter.AdjacencyTransformation import ATPureAdjacency
import argparse

ADJACENCIES = {
    'pure': ATPureAdjacency,
}

def main():
    parser = argparse.ArgumentParser(
        prog='PyReporter',
        description='PyReporter tool for metabolomics predictions',
    )
    parser.add_argument('expression_file')

    parser.add_argument('-b', '--baseline')
    parser.add_argument('-m', '--model')
    parser.add_argument('-a', '--adjacency')
    parser.add_argument('-e', '--count')
    parser.add_argument('-v', '--verbose')

    args = parser.parse_args()

    adjacency = ADJACENCIES.get(args.adjacency, ATPureAdjacency)
    
    
    # results = workflow_Fang2012(
    #     cobra_model,
    #     mapped_genes_baseline,
    #     mapped_genes_comparison,
    #     adjacency = ADJACENCIES.get(args.adjacency),
    #     ranking = pr.PageRank.PageRankNX,
    #     gene_fill = 1.,
    # ) -> pd.Series: