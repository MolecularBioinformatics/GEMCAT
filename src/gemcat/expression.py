"""
Module containing expression data integration
"""

import abc
import logging
import re
from typing import Optional, Tuple

import cobra
import numpy as np
import pandas as pd
from sympy.parsing.sympy_parser import parse_expr

from . import utils, verification

geomean = utils.geometric_mean


def read_simple_gpr_from_cobra(
    model: cobra.Model,
) -> dict[str, list[str]]:
    """
    Reads a simple version of GPR from a COBRA model:
    All reactions are mapped to all involved genes.
    :return: dictionary with keys reactions, values associated genes
    :rtype: dict[str: list[str]]
    """
    genes = []
    for rxn in model.reactions:
        genes.append([g.id for g in rxn.genes])

    genes_dict = dict(zip([rxn.id for rxn in model.reactions], genes))
    return genes_dict


def read_gpr_strings_from_cobra(
    model: cobra.Model,
) -> Tuple[dict[str, str], list[str]]:
    """
    Extracts the gene product rules (GPRs) from a model.
    :param model: Model from which to extract GPR
    :type model: cobra.Model
    :return: Tuple of extracted GPR (dict) and list of all genes.
    :rtype: Tuple[ dict[str, str], list[str] ]
    """
    gpr = {rxn.id: rxn.gene_reaction_rule for rxn in model.reactions}
    genes = {rxn.id: [gene.id for gene in rxn.genes] for rxn in model.reactions}
    return gpr, genes


def fillna_mean(series: pd.Series) -> pd.Series:
    """
    Fills all NaN values in self.mapped_values with the mean
    of the other data.
    """
    return series.fillna(series.mean())


class ExpressionIntegration(abc.ABC):
    """
    Abstract base class for Expression to define interface.
    Do NOT use.
    """

    @abc.abstractmethod
    def __init__(self, data: pd.Series, gpr: dict[str, list[str]], gene_fill: float):
        self.mapped_values: Optional[pd.Series] = None

    @staticmethod
    def _verify_data(data):
        """
        Verify expression data
        """
        verification.raise_for_duplicated_index(data)
        verification.raise_for_non_int_float_64_dtype(data)
        data = verification.convert_index_to_str(data)
        return data

    @abc.abstractmethod
    def load_gpr(self, gpr: dict[str, list[str]]):
        """
        Load GPR into usable format
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def map(self):
        """
        Map values onto GPR
        """
        raise NotImplementedError()

    def get_mapped_values(self) -> np.array:
        """
        Returns the mapped values as a NumPy array.
        :return: Gene expression values
        :rtype: 1-D np.array (r)
        """
        return self.mapped_values.values

    def fillna(self, fill_fn=lambda series: series.fillna(1.0)) -> None:
        """
        Fills NaN values in the mapped reaction values according to a function passed.
        The passed function receives a pd.Series object and returns a filled version.
        :param fill_fn: Function filling NaN values in a pd.Series
        :type fill_fn: Function (fill_fn:: pd.Series -> pd.Series)
        """
        self.mapped_values = fill_fn(self.mapped_values)


class GeometricAndAverageMeans(ExpressionIntegration):
    """
    Integration of gene expression based on Fang et al. (2012)
    doi.org/10.1371/journal.pcbi.1002688

    Where reactions depend on multiple genes being active (logical AND),
    the geometric mean of the values for each gene is used.
    Where reactions depend on either of multiple genes (logical OR),
    the gene values are added up.

    """

    def __init__(
        self,
        gpr: dict[str, str],
        reaction_gene_mapping: dict[str, list[str]],
        data: pd.Series,
        gene_fill=0.0,
    ):
        """
        Create a component for the integration of expression data.
        This adheres to the algorithm laid out by Fang et al. 2012.
        :param gpr: Gene product rule for each reaction, dict[Reaction: GPRstring]
        :type gpr: dict[str, str]
        :param reaction_gene_mapping: Mapping which genes participate
        in which reactions, dict[reaction: list[geneString]]
        :type reaction_gene_mapping: dict[str, list[str]]
        :param data: Omics data mapped to genes,
        pd.Series[geneString: value]
        :type data: pd.Series
        :param gene_fill: value to fill in for genes missing in the data set
        :type gene_fill: float
        """
        self.quant_gpr: pd.Series = None
        self.mapped_values: pd.Series = None
        self.gene_fill = gene_fill
        self.data = self._verify_data(data)
        self.load_gpr(gpr)
        self.reaction_gene_list = reaction_gene_mapping
        self.rewrite_gpr()
        self.map()
        self.fillna()

    def load_gpr(self, gpr: dict[str, str]):
        """
        Load GPRs into an appropriate format
        :param gpr: GPR
        :type gpr: dict[str, str]
        """
        self.gpr = pd.Series(gpr)

    def rewrite_gpr(self):
        """
        Rewrite all GPRs to evaluable strings
        """
        quant_gpr = self.gpr.index.map(self.rewrite_single_gpr)
        self.quant_gpr = pd.Series(list(quant_gpr), self.gpr.index)

    def rewrite_single_gpr(self, rxn: str) -> str:
        """
        Rewrite a single GPR to an evaluable string that describes
        the calculation of its value.
        :param rxn: reaction ID
        :type rxn: str
        :return: Rewritten GPR
        :rtype: str
        """
        gpr = self.gpr[rxn]
        genes = self.reaction_gene_list.get(rxn, "")
        if not gpr:
            return ""
        if not isinstance(gpr, str):
            return ""

        for gene in genes:
            gene_str = str(gene)
            gene_val = float(self.data.get(gene_str, self.gene_fill))
            gpr = gpr.replace(gene_str, f"{gene_val}")

        gpr = gpr.replace("or", "+")
        re_float = "\d*\.\d*"
        re_and = re.compile(f"{re_float}(?: and {re_float})+")
        hits = re_and.findall(gpr)
        if not hits:
            return gpr
        for hit in hits:
            gids = hit.split("and")
            gids = [x.strip() for x in gids]
            gids = ", ".join(gids)
            new = f"geomean({gids})"
            gpr = gpr.replace("(" + hit + ")", "(" + new + ")")

        return gpr

    def map(self):
        self.mapped_values = self.quant_gpr.map(self.eval_single_gpr)

    @staticmethod
    def eval_single_gpr(
        gpr: str,
    ):
        """
        Evaluate a gene product rule (GPR)
        :param gpr: GPR string
        :type gpr: str
        :return: _description_
        :rtype: _type_
        """
        if len(gpr) == 0:
            return np.nan
        result = parse_expr(gpr, {"geomean": geomean})
        try:
            return float(result)
        except TypeError:
            logging.debug("Failed: %s", gpr)
            return np.nan


class ExpressionMapSingleAverage(ExpressionIntegration):
    """
    Simple integration of expression data by mapping each reaction
    onto the mean expression of the genes involved in it.
    """

    def __init__(
        self,
        data: pd.Series,
        gpr: dict[str, list[str]],
        gene_fill=0.0,
    ):
        """
        Initialize object and map expression data onto reactions.
        :param data: Gene expression data
        :type data: pd.Series
        :param gpr: Gene product rules
        :type gpr: dict[str, list[str]]
        """
        self.gpr = None
        self.load_gpr(gpr)
        self.mapped_values = None
        data = self._verify_data(data)
        self.map(data)

    def load_gpr(
        self,
        gpr: dict[str, list[str]],
    ):
        """
        Load gene product rule into the object.
        :param gpr: [description]
        :type gpr: dict[str, list[str]]
        """
        self.gpr = gpr

    def map(self, data: pd.Series):
        """
        Map the expression data using the gene product rule.
        Save mapped values in object as a pandas series.
        :param data: Gene expression data
        :type data: pd.Series
        :raises ValueError: Raised if no GPR is set on the object
        """
        if not self.gpr:
            err = "Attempting to map expression data while GPR has not yet been set."
            logging.error(err)
            raise ValueError(err)
        vals_dict = {}
        for reaction, genes in self.gpr.items():
            expression = []
            for gene in genes:
                try:
                    expression.append(data.loc[gene])
                except KeyError:
                    expression.append(0)
            if not expression:
                vals_dict[reaction] = np.nan
            else:
                vals_dict[reaction] = sum(expression) / len(expression)
        self.mapped_values = pd.Series(vals_dict)
        self.fillna()
