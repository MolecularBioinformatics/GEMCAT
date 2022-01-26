#!/usr/bin/python

import pandas as pd
import cobra
import pyreporter as pr


regex_recon3d_old = pr.regexes.REGEX['Recon3D_old']

def workflow_single(
    cobra_model: cobra.Model, 
    mapped_genes: pd.Series,
    adjacency = pr.AdjacencyTransformation.ATPureAdjacency,
    ranking = pr.PageRank.PageRankNX,
    ) -> pd.Series:
    
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
    re_gene = regex_recon3d_old,
    adjacency = pr.AdjacencyTransformation.ATPureAdjacency,
    ranking = pr.PageRank.PageRankNX,
    ) -> pd.Series:

    model = pr.io.convert_cobra_model(cobra_model)
    model.AT = adjacency
    model.ranking = ranking
    ex_baseline = pr.Expression.ExpressionFang2012(
        model = cobra_model,
        data = mapped_genes_comparison,
        re_gene = re_gene
    )
    ex_comparison = pr.Expression.ExpressionFang2012(
        model = cobra_model,
        data = mapped_genes_baseline,
        re_gene = re_gene
    )
    model.load_expression(ex_comparison)
    results_comparison = model.calculate()
    model.load_expression(ex_baseline)
    results_baseline = model.calculate()
    results = results_comparison / results_baseline

    return results