#!/usr/bin/python

import pandas as pd
import cobra
import pyreporter as pr


def workflow_single(
    cobra_model: cobra.Model, 
    mapped_genes: pd.Series,
    adjacency = pr.AdjacencyTransformation.ATPureAdjacency,
    ranking = pr.PageRank.PageRankNX,
    ) -> pd.Series:
    """
    Workflow using Expression data integration via simple averaging, provides only the output from a single expression data set.
    :param cobra_model: cobra Model used as a base for modeling
    :type cobra_model: cobra.Model
    :param mapped_genes: Levels of gene/protein expression values
    :type mapped_genes: pd.Series
    :param adjacency: Adjacency metric, defaults to pr.AdjacencyTransformation.ATPureAdjacency
    :type adjacency: pr.AdjacencyTransformation.AdjacencyTransformation, optional
    :param ranking: Ranking algorithm class, defaults to pr.PageRank.PageRankNX
    :type ranking: pr.PageRank.Ranking, optional
    :return: Normalized relative metabolite scores: comparison / baseline
    :rtype: pd.Series
    """
    model = pr.io.convert_cobra_model(cobra_model)
    model.AT = adjacency
    model.ranking = ranking
    gpr = pr.Expression.read_simple_gpr_from_cobra(cobra_model)
    ex = pr.Expression.ExpressionMapSingleAverage(mapped_genes, gpr)
    model.load_expression(ex)
    results = model.calculate()
    return results, model

def workflow_ratio(
    cobra_model: cobra.Model,
    mapped_genes_baseline: pd.Series,
    mapped_genes_comparison: pd.Series,
    adjacency = pr.AdjacencyTransformation.ATPureAdjacency,
    ranking = pr.PageRank.PageRankNX,
    ) -> pd.Series:
    """
    Workflow for differerential expression between two datasets.
    Integrates expression data via simple averaging.
    :param cobra_model: cobra Model used as a base for modeling
    :type cobra_model: cobra.Model
    :param mapped_genes_baseline: Baseline levels of gene/protein expression values
    :type mapped_genes_baseline: pd.Series
    :param mapped_genes_comparison: Comparison levels of gene/protein expression values
    :type mapped_genes_comparison: pd.Series
    :param adjacency: Adjacency metric, defaults to pr.AdjacencyTransformation.ATPureAdjacency
    :type adjacency: pr.AdjacencyTransformation.AdjacencyTransformation, optional
    :param ranking: Ranking algorithm class, defaults to pr.PageRank.PageRankNX
    :type ranking: pr.PageRank.Ranking, optional
    :return: Normalized relative metabolite scores: comparison / baseline
    :rtype: pd.Series
    """
    model = pr.io.convert_cobra_model(cobra_model)
    model.AT = adjacency
    model.ranking = ranking
    gpr = pr.Expression.read_simple_gpr_from_cobra(cobra_model)
    ex = pr.Expression.ExpressionMapSingleAverage(mapped_genes_comparison, gpr)
    ex_baseline = pr.Expression.ExpressionMapSingleAverage(mapped_genes_baseline, gpr)
    model.load_expression(ex)
    results = model.calculate()
    model.load_expression(ex_baseline)
    results_baseline = model.calculate()
    results = results / results_baseline

    return results

def workflow_Fang2012(
    cobra_model: cobra.Model,
    mapped_genes_baseline: pd.Series,
    mapped_genes_comparison: pd.Series,
    adjacency = pr.AdjacencyTransformation.ATPureAdjacency,
    ranking = pr.PageRank.PageRankNX,
    gene_fill = 0.,
    ) -> pd.Series:
    """
    Standard workflow integrating expression data via the method provided by Fang2012.
    :param cobra_model: COBRA model to be used to generate PyReporter model and GPR
    :type cobra_model: cobra.Model
    :param mapped_genes_baseline: Baseline levels of gene/protein expression values
    :type mapped_genes_baseline: pd.Series
    :param mapped_genes_comparison: Comparison levels of gene/protein expression values
    :type mapped_genes_comparison: pd.Series
    :param adjacency: Adjacency metric, defaults to pr.AdjacencyTransformation.ATPureAdjacency
    :type adjacency: pr.AdjacencyTransformation.AdjacencyTransformation, optional
    :param ranking: Ranking algorithm class, defaults to pr.PageRank.PageRankNX
    :type ranking: pr.PageRank.Ranking, optional
    :return: Normalized relative metabolite scores: comparison / baseline
    :rtype: pd.Series
    """
    model = pr.io.convert_cobra_model(cobra_model)
    gpr, rxn_gene_mapping = pr.Expression.read_gpr_strings_from_cobra(cobra_model)
    model.AT = adjacency
    model.ranking = ranking
    ex_baseline = pr.Expression.ExpressionFang2012(
        gpr = gpr,
        reaction_gene_mapping = rxn_gene_mapping,
        data = mapped_genes_baseline,
        gene_fill = gene_fill,
    )
    ex_comparison = pr.Expression.ExpressionFang2012(
        gpr = gpr,
        reaction_gene_mapping = rxn_gene_mapping,
        data = mapped_genes_comparison,
        gene_fill = gene_fill,
    )
    model.load_expression(ex_comparison)
    results_comparison = model.calculate()
    model.load_expression(ex_baseline)
    results_baseline = model.calculate()
    results = results_comparison / results_baseline

    return (results, results_comparison, results_baseline)
