#!/usr/bin/python

"""
Module containing complete end-to-end workflows in Python.
Usable as documentation, or as quick runs from within a Python script.
"""

import cobra
import pandas as pd

from . import adjacency_transformation as at
from . import expression as ex
from . import io
from . import ranking as pr


def workflow_avg_single(
    cobra_model: cobra.Model,
    mapped_genes: pd.Series,
) -> pd.Series:
    """
    Workflow using Expression data integration via simple averaging,
    provides only the output from a single expression data set.
    :param cobra_model: cobra Model used as a base for modeling
    :type cobra_model: cobra.Model
    :param mapped_genes: Levels of gene/protein expression values
    :type mapped_genes: pd.Series
    :param adjacency: Adjacency metric, defaults to ATPureAdjacency
    :type adjacency: at.adjacency_transformation, optional
    :param ranking: Ranking algorithm class, defaults to PagerankNX
    :type ranking: pr.Ranking, optional
    :return: Normalized metabolite scores
    :rtype: pd.Series
    """
    model = io.convert_cobra_model(cobra_model)
    model.adjacency_transformation = at.ATPureAdjacency()
    model.ranking = pr.PagerankNX()
    gpr = ex.read_simple_gpr_from_cobra(cobra_model)
    expression = ex.ExpressionMapSingleAverage(mapped_genes, gpr)
    model.load_expression(expression)
    results = model.calculate()
    return results


def workflow_avg_ratio(
    cobra_model: cobra.Model,
    mapped_genes_baseline: pd.Series,
    mapped_genes_comparison: pd.Series,
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
    :param adjacency: Adjacency metric, defaults to at.ATPureAdjacency
    :type adjacency: at.adjacency_transformation, optional
    :param ranking: Ranking algorithm class, defaults to pr.PagerankNX
    :type ranking: pr.Ranking, optional
    :return: Normalized relative metabolite scores: comparison / baseline
    :rtype: pd.Series
    """
    model = io.convert_cobra_model(cobra_model)
    model.adjacency_transformation = at.ATPureAdjacency()
    model.ranking = pr.PagerankNX()
    gpr = ex.read_simple_gpr_from_cobra(cobra_model)
    expression = ex.ExpressionMapSingleAverage(mapped_genes_comparison, gpr)
    expression_baseline = ex.ExpressionMapSingleAverage(mapped_genes_baseline, gpr)
    model.load_expression(expression)
    results = model.calculate()
    model.load_expression(expression_baseline)
    results_baseline = model.calculate()
    results = results / results_baseline

    return results


def workflow_standard(
    cobra_model: cobra.Model,
    mapped_genes_baseline: pd.Series,
    mapped_genes_comparison: pd.Series,
    gene_fill=1.0,
) -> pd.Series:
    """
    Standard workflow integrating expression data via the method provided by Fang2012.
    :param cobra_model: COBRA model to be used to generate PyReporter model and GPR
    :type cobra_model: cobra.Model
    :param mapped_genes_baseline: Baseline levels of gene/protein expression values
    :type mapped_genes_baseline: pd.Series
    :param mapped_genes_comparison: Comparison levels of gene/protein expression values
    :type mapped_genes_comparison: pd.Series
    :param adjacency: Adjacency metric, defaults to at.ATPureAdjacency
    :type adjacency: at.adjacency_transformation, optional
    :param ranking: Ranking algorithm class, defaults to pr.PagerankNX
    :type ranking: pr.Ranking, optional
    :param gene_fill: Value to fill in for genes missing in input data, defaults to 0.
    :type gene_fill: float
    :return: Normalized relative metabolite scores: comparison / baseline
    :rtype: pd.Series
    """
    model = io.convert_cobra_model(cobra_model)
    gpr, rxn_gene_mapping = ex.read_gpr_strings_from_cobra(cobra_model)
    model.adjacency_transformation = at.ATPureAdjacency()
    model.ranking = pr.PagerankNX()
    ex_baseline = ex.GeometricAndAverageMeans(
        gpr=gpr,
        reaction_gene_mapping=rxn_gene_mapping,
        data=mapped_genes_baseline,
        gene_fill=gene_fill,
    )
    ex_comparison = ex.GeometricAndAverageMeans(
        gpr=gpr,
        reaction_gene_mapping=rxn_gene_mapping,
        data=mapped_genes_comparison,
        gene_fill=gene_fill,
    )
    model.load_expression(ex_comparison)
    results_comparison = model.calculate()
    model.load_expression(ex_baseline)
    results_baseline = model.calculate()
    results = results_comparison / results_baseline

    return results
