from warnings import warn
from typing import List, Optional
import numpy as np
import pyreporter.PageRank as PR
import pyreporter.AdjacencyTransformation as AT
import pyreporter.utils as utils
import pandas as pd
import pyreporter.Expression as EX
import networkx as nx

class Model():
    """
    Central class of the analysis framework. 
    Collects essential features and algorithmic interfaces.
    Implements the central work flow.
    
    Public attributes:
    S - stoichiometric matrix (m x r)
    ranking - the ranking algorithm used (currently only PageRank)
    metabolite_names - list of metabolite names in order
    expression - vector of expression values (1 x r)
    scores - metabolite scores used for personalization

    Public methods:
    calculate - calculates metabolite scores with current values
    load expression - load expression data
    """

    def __init__(
        self, 
        S: np.array,
        metabolite_names: List[str],
        reversibilities: List[bool],
        at: Optional[AT.ATPureAdjacency] = None,
        ranking: Optional[PR.PageRankNX] = None,
        metabolite_seeds: Optional[List[float]] = None
        ):
        """
        Create a model object
        :param S: stoichiometric matrix, including reversible reactions
        :type S: np.array
        :param metabolite_names: List of metabolite names in order
        :type metabolite_names: List[str]
        :param reversibilities: List of reaction reversibilities
        :type reversibilities: List[bools]
        """
        self.S = S
        self.dimensions = self.S.shape
        self.expression_shape = (1, self.dimensions[1])
        
        if at is None:
            at = AT.ATPureAdjacency()
        self.AT = at
        if ranking is None:
            ranking = PR.PageRankNX()
        self.ranking = ranking

        self.metabolite_names = metabolite_names
        self.expression_vector = None
        self.reversibilities = reversibilities
        self.expression = None
        self._update_expression_vector()
        self._A_is_current = False
        self.scores = None
        self.seeds = None
        self.load_metabolite_seeds(metabolite_seeds)

    def load_metabolite_seeds(self, seeds):
        """
        Load metabolite seeds and convert to numpy array if necessary.
        :param seeds: Metabolite score seeds
        :type seeds: List of Numpy array
        :raises ValueError: In case of incompatible dimensions
        """
        if seeds is None:
            self.seeds = None
            return
        if not isinstance(seeds, list):
            raise TypeError('Expected metabolite seeds to be of type list')
        if not (len(seeds) == self.dimensions[0]):
            raise ValueError("Length of seeds must be equal to number of metabolites")
        self.seeds = seeds

    def _update_A(self):
        """
        Calculate the adjacency matrix with currently set values 
        and store it in the model.
        :return: Adjacency matrix
        :rtype: np.array (m x m)
        """
        self.A = self.AT.transform(
            self.S, 
            self.reversibilities, 
            self.expression_vector
        )

    def load_expression(self, expression: EX.Expression):
        """
        Load expression data into the model.
        Expression data needs to be in order matching S.
        :param omics_array: Array of reaction scores 
        :type omics_array: np.array (r)
        """
        if not isinstance(expression, EX.Expression):
            raise TypeError('Needs to be an Expression object')
        if self.expression:
            warn('Previous expression data overwritten')
        self._A_is_current = False
        self.expression = expression
        self._update_expression_vector()

    def calculate(
        self, 
        graph_args: Optional[dict] = None, 
        pr_args: Optional[dict] = None
        ) -> pd.Series:
        """
        Calculate scores with current S, expression, and metabolite score seeds.
        :param graph_args: Arguments to pass to graph creation, defaults to {}
        :type graph_args: dict, optional
        :param pr_args: Args to pass to ranking, defaults to {}
        :type pr_args: dict, optional
        :return: Scores for each metabolite
        :rtype: pd.Series
        """
        if graph_args is None:
            graph_args = {}
        if pr_args is None:
            pr_args = {}
        self._check_and_reload__A()
        scores = self.ranking.propagate(
            self.A, 
            self.seeds,
            self.metabolite_names,
            graph_args, 
            pr_args
        )
        self.scores = scores
        return utils._annotate(scores, self.metabolite_names)

    def _update_expression_vector(self):
        """
        Initializes expression scores
        """
        if self.expression:
            self.expression_vector = self.expression.get_mapped_values()
        else:
            self.expression_vector = np.ones(
                shape = self.expression_shape,
                dtype = np.float64,
                )
    
    def _check_and_reshape_expression_vector(self):
        """
        Checks shape of expression vector. 
        If it doesn't fit S, reshapes it to fit.
        :raises ValueError: Raised if the expression vector is the wrong shape
        """
        if not self.expression_vector.shape == self.expression_shape:
            try:
                self.expression_vector = self.expression_vector.reshape(self.expression_shape)
            except:
                raise ValueError('Current expression vector is wrong length')
    
    def _check_and_reload__A(self):
        """
        Checks whether A is current, if not, reloads it.
        """
        if not self._A_is_current:
            self._update_A()
            self._A_is_current = True

    def get_subnetworks(self):
        """
        Returns subnetworks (weakly connected) in the current A.
        :return: 
        :rtype: [type]
        """
        self._check_and_reload__A()
        g = nx.DiGraph(self.A)
        wcc = nx.algorithms.weakly_connected_components(g)
        wcc_full = []
        for subnetwork in wcc:
            sub_mets = []
            for met_no in subnetwork:
                sub_mets.append(self.metabolite_names[met_no])
            wcc_full.append(sub_mets)
        return wcc_full