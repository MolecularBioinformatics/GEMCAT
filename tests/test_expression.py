import numpy as np
import pandas as pd
from fixtures import models

from gemcat import expression as ex

R_TOLERANCE = 10**-5
gpr_mini = {
    "R1": ["G1"],
    "R2": ["G2"],
    "R3": ["G3"],
    "R4": ["G4"],
}
genes = [f"G{i}" for i in range(1, 9)]
genes_gid = [f"000{i}.{i}" for i in range(1, 9)]
gene_vals = [float(i) for i in range(1, 9)]
genes_mini_complex = pd.Series(dict(zip(genes, gene_vals)))
gids_mini_complex = pd.Series(dict(zip(genes_gid, gene_vals)))


def test_read_simple_gpr_from_cobra(models):
    gpr = ex.read_simple_gpr_from_cobra(models["mini"])
    assert gpr == gpr_mini


def test_ExpressionMapSingleAverage(models):
    gene_map = pd.Series({"G1": 2.0, "G2": 2.0, "G3": 2.0, "G4": 2.0})
    gpr = ex.read_simple_gpr_from_cobra(models["mini"])
    exp = ex.ExpressionMapSingleAverage(gene_map, gpr)
    expected = np.array([2, 2, 2, 2])
    assert np.allclose(exp.get_mapped_values(), expected)


def test_ExpressionMapFang2012_simple(models):
    gpr, gene_dict = ex.read_gpr_strings_from_cobra(models["mini_complex_gpr"])
    exp = ex.GeometricAndAverageMeans(
        gpr=gpr,
        reaction_gene_mapping=gene_dict,
        data=genes_mini_complex,
    )
    assert isinstance(exp.gpr, pd.Series)
    assert np.isclose(exp.mapped_values["R4"], 4.0, rtol=R_TOLERANCE)


def test_ExpressionMapFang2012_complex(models):
    gpr, gene_dict = ex.read_gpr_strings_from_cobra(models["mini_complex_gpr"])
    exp = ex.GeometricAndAverageMeans(
        gpr=gpr,
        reaction_gene_mapping=gene_dict,
        data=genes_mini_complex,
    )
    assert isinstance(exp.gpr, pd.Series)
    assert np.isclose(exp.mapped_values["R3"], 15.4157, rtol=R_TOLERANCE)


def test_ExpressionMapFang2012_gids(models):
    gpr, gene_dict = ex.read_gpr_strings_from_cobra(models["mini_complex_gpr_gids"])
    exp = ex.GeometricAndAverageMeans(
        gpr,
        gene_dict,
        gids_mini_complex,
    )
    assert isinstance(exp.gpr, pd.Series)
    assert np.isclose(exp.mapped_values["R3"], 15.4157, rtol=R_TOLERANCE)


def test_read_complex_gpr_structure(models):
    gpr = ex.read_gpr_strings_from_cobra(models["mini_complex_gpr"])
    assert isinstance(gpr, tuple)
    assert isinstance(gpr[0], dict)
    assert isinstance(gpr[1], dict)


def test_read_complex_gpr_gene_rxn(models):
    rxn_genes = ex.read_gpr_strings_from_cobra(models["mini_complex_gpr"])[1]
    expected = ["G3", "G4", "G5", "G6", "G7", "G8"]
    assert set(rxn_genes["R3"]) == set(expected)


def test_read_complex_gpr_gpr(models):
    expected = models["mini_complex_gpr"].reactions.get_by_id("R3").gene_reaction_rule
    gpr = ex.read_gpr_strings_from_cobra(models["mini_complex_gpr"])[0]
    assert gpr["R3"] == expected
