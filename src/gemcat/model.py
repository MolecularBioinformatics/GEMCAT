#!/usr/bin/python

"""
Model structure central to the framework
"""

import logging
from typing import Optional

import networkx as nx
import numpy as np
import pandas as pd

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
    S - stoichiometric matrix (m x r)
    ranking - the ranking algorithm used (currently only Pagerank)
    metabolite_names - list of metabolite names in order
    expression - vector of expression values (1 x r)
    scores - metabolite scores used for personalization

    Public methods:
    calculate - calculates metabolite scores with current values
    load expression - load expression data
    """

    def __init__(
        self,
        stoichiometric_matrix: np.array,
        metabolite_names: list[str],
        reversibilities: list[bool],
        adjacency: Optional[at.AdjacencyTransformation] = None,
        ranking: Optional[pr.Ranking] = None,
        metabolite_seeds: Optional[list[float]] = None,
    ):
        """
        Create a model object
        :param stoichiometric_matrix: stoichiometric matrix,
        including reversible reactions
        :type stoichiometric_matrix: np.array [m x r]
        :param metabolite_names: List of metabolite names in order
        :type metabolite_names: list[str]
        :param reversibilities: List of reaction reversibilities
        :type reversibilities: list[bool]
        :param adjacency: Adjacency calculation object, defaults to using PureAdjacency
        :type adjacency: Optional[AT.ATPureAdjacency], optional
        :param ranking: _description_, defaults to using PageRankNX
        :type ranking: Optional[PR.PagerankNX], optional
        :param metabolite_seeds: _description_, defaults to None
        :type metabolite_seeds: Optional[list[float]], optional
        """
        self.stoichiometric_matrix = stoichiometric_matrix
        self.adjacencies: np.ndarray = None
        self.dimensions = self.stoichiometric_matrix.shape
        self.expression_shape = (1, self.dimensions[1])

        if adjacency is None:
            adjacency = at.ATPureAdjacency()
        self.adjacency_transformation = adjacency
        if ranking is None:
            ranking = pr.PagerankNX()
        self.ranking = ranking

        self.metabolite_names = metabolite_names
        self.expression_vector = None
        self.reversibilities = reversibilities
        self.expression = None
        self._update_expression_vector()
        self._adjacencies_are_current = False
        self.scores = None
        self.seeds = None
        self.load_metabolite_seeds(metabolite_seeds)

    def load_metabolite_seeds(self, seeds):
        """
        Load metabolite seeds and convert to numpy array if necessary.
        :param seeds: Metabolite score seeds
        :type seeds: list[np.ndarray]
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

    def _update_adjacencies(self):
        """
        Calculate the adjacency matrix with currently set values
        and store it in the model.
        :return: Adjacency matrix
        :rtype: np.array (m x m)
        """
        self.adjacencies = self.adjacency_transformation.transform(
            self.stoichiometric_matrix, self.reversibilities, self.expression_vector
        )

    def load_expression(self, expression: ex.ExpressionIntegration):
        """
        Load expression data into the model.
        Expression data needs to be in order matching S.
        :param omics_array: Array of reaction scores
        :type omics_array: np.array (r)
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
        :type graph_args: dict, optional
        :param pr_args: Args to pass to ranking, defaults to None into {}
        :type pr_args: dict, optional
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

    def _update_expression_vector(self):
        """
        Initializes expression scores
        """
        if self.expression:
            self.expression_vector = self.expression.get_mapped_values()
        else:
            self.expression_vector = np.ones(
                shape=self.expression_shape,
                dtype=np.float64,
            )

    def _check_and_reshape_expression_vector(self):
        """
        Checks shape of expression vector.
        If it doesn't fit S, reshapes it to fit.
        :raises ValueError: Raised if the expression vector is the wrong shape
        """
        if self.expression_vector.shape == self.expression_shape:
            return
        try:
            self.expression_vector = self.expression_vector.reshape(
                self.expression_shape
            )
        except Exception as original_err:
            err = "Current expression vector is wrong length"
            logging.error(err)
            raise ValueError(err) from original_err

    def _check_and_reload_adjacencies(self):
        """
        Checks whether A is current, if not, reloads it.
        """
        if not self._adjacencies_are_current:
            self._update_adjacencies()
            self._adjacencies_are_current = True

    def get_subnetworks(self):
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
