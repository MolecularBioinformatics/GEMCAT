#!/usr/bin/python

"""
Model structure central to the framework
"""

import logging
from typing import Optional

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import adjacency_transformation as at
from . import expression as ex
from . import ranking as pr
from . import utils


class Model:
    """
    Central class of the analysis framework.
    Collects essential features and algorithmic interfaces.
    Implements the central work flow.

    Public attributes:
    S - sparse stoichiometric matrix (m x r)
    ranking - the ranking algorithm used (currently only Pagerank)
    metabolite_names - list of metabolite names in order
    expression - vector of expression values (r,)
    scores - metabolite scores used for personalization

    Public methods:
    calculate - calculates metabolite scores with current values
    load expression - load expression data
    """

    def __init__(
        self,
        stoichiometric_matrix: sp.sparray,
        metabolite_names: list[str],
        reversibilities: list[bool],
        adjacency: Optional[at.AdjacencyTransformation] = None,
        ranking: Optional[pr.Ranking] = None,
        metabolite_seeds: Optional[list[float]] = None,
    ):
        """
        Create a model object
        :param stoichiometric_matrix: sparse stoichiometric matrix,
        including reversible reactions (m, r)
        :type stoichiometric_matrix: sp.sparray
        :param metabolite_names: List of metabolite names in order
        :type metabolite_names: list[str]
        :param reversibilities: List of reaction reversibilities
        :type reversibilities: list[bool]
        :param adjacency: Adjacency calculation object, an instance rather than a
        class, defaults to ATPureAdjacency
        :type adjacency: Optional[at.AdjacencyTransformation], optional
        :param ranking: Ranking algorithm object, defaults to PagerankNX
        :type ranking: Optional[pr.Ranking], optional
        :param metabolite_seeds: Per-metabolite personalization weights,
        defaults to None
        :type metabolite_seeds: Optional[list[float]], optional
        """
        utils.require_sparse(stoichiometric_matrix, "the stoichiometric matrix")
        self.stoichiometric_matrix: sp.sparray = stoichiometric_matrix
        self.adjacencies: Optional[sp.csr_array] = None
        self.dimensions = self.stoichiometric_matrix.shape
        self.expression_shape = (self.dimensions[1],)

        if adjacency is None:
            adjacency = at.ATPureAdjacency()
        self.adjacency_transformation = adjacency
        if ranking is None:
            ranking = pr.PagerankNX()
        self.ranking = ranking

        self.metabolite_names = metabolite_names
        self.reversibilities = reversibilities
        self.expression: Optional[ex.ExpressionIntegration] = None
        self.expression_vector: np.ndarray
        self._update_expression_vector()
        self._adjacencies_are_current = False
        self.scores: Optional[np.ndarray] = None
        self.seeds: Optional[list[float]] = None
        self.load_metabolite_seeds(metabolite_seeds)

    def load_metabolite_seeds(self, seeds: Optional[list[float]]) -> None:
        """
        Load metabolite seeds to use as PageRank personalization weights.
        :param seeds: Metabolite score seeds, one per metabolite
        :type seeds: Optional[list[float]]
        :raises TypeError: If seeds is neither None nor a list
        :raises ValueError: In case of incompatible dimensions
        """
        if seeds is None:
            self.seeds = None
            return
        if not isinstance(seeds, list):
            raise TypeError("Expected metabolite seeds to be of type list")
        if len(seeds) != self.dimensions[0]:
            raise ValueError("Length of seeds must be equal to number of metabolites")
        self.seeds = seeds

    def _update_adjacencies(self) -> None:
        """
        Calculate the adjacency matrix from the current values and store it in
        self.adjacencies.
        """
        self.adjacencies = self.adjacency_transformation.transform(
            self.stoichiometric_matrix, self.reversibilities, self.expression_vector
        )

    def load_expression(self, expression: ex.ExpressionIntegration) -> None:
        """
        Load expression data into the model.
        Expression data needs to be in order matching S.
        :param expression: Expression integration; one score per reaction
        :type expression: ex.ExpressionIntegration
        :raises TypeError: If the argument is not an ExpressionIntegration
        """
        if not isinstance(expression, ex.ExpressionIntegration):
            raise TypeError("Needs to be an Expression object")
        if self.expression:
            logging.debug("Previous expression data overwritten")
        self._adjacencies_are_current = False
        self.expression = expression
        self._update_expression_vector()

    def calculate(
        self, graph_args: Optional[dict] = None, pr_args: Optional[dict] = None
    ) -> pd.Series:
        """
        Calculate scores with current S, expression, and metabolite score seeds.
        :param graph_args: Arguments to pass to graph creation, defaults to None into {}
        :type graph_args: Optional[dict]
        :param pr_args: Args to pass to ranking, defaults to None into {}
        :type pr_args: Optional[dict]
        :return: Scores for each metabolite
        :rtype: pd.Series
        """
        if graph_args is None:
            graph_args = {}
        if pr_args is None:
            pr_args = {}
        self._check_and_reload_adjacencies()
        scores = self.ranking.propagate(
            self.adjacencies, self.seeds, self.metabolite_names, graph_args, pr_args
        )
        self.scores = scores
        return utils.annotate_scores(scores, self.metabolite_names)

    def _update_expression_vector(self) -> None:
        """
        Rebuild the per-reaction expression vector from the loaded expression data.

        With no expression loaded it falls back to all ones. Either way
        the vector is checked and flattened to (r,) for the transforms.
        """
        if self.expression:
            raw = self.expression.get_mapped_values()
        else:
            raw = np.ones(shape=self.expression_shape, dtype=np.float64)
        self.expression_vector = utils.as_reaction_vector(raw, self.dimensions[1])

    def _check_and_reshape_expression_vector(self) -> None:
        """
        Normalizes the current expression vector to the shape S expects.
        :raises ValueError: Raised if the expression vector is the wrong length
        """
        self.expression_vector = utils.as_reaction_vector(
            self.expression_vector, self.dimensions[1]
        )

    def _check_and_reload_adjacencies(self) -> None:
        """
        Checks whether A is current, if not, reloads it.
        """
        if not self._adjacencies_are_current:
            self._update_adjacencies()
            self._adjacencies_are_current = True

    def get_subnetworks(self) -> list[list[str]]:
        """
        Returns subnetworks (weakly connected) in the current A.
        :return: List of weakly connected subnetworks
        (names of metabolites in the subnet)
        :rtype: list[list[str]]
        """
        self._check_and_reload_adjacencies()
        graph = nx.DiGraph(self.adjacencies)
        connected = nx.algorithms.weakly_connected_components(graph)
        connected_full = []
        for subnetwork in connected:
            sub_mets = [self.metabolite_names[met_no] for met_no in subnetwork]
            connected_full.append(sub_mets)
        return connected_full
