#!/usr/bin/python

"""
Algorithms to propagate scores from the adjacency matrix
"""

import logging
from abc import ABC, abstractmethod
from typing import Dict, Optional

import networkx as nx
import numpy as np


class Ranking(ABC):
    """
    Abstract base class for Ranking methods. Do NOT use.
    """

    @staticmethod
    @abstractmethod
    def propagate(
        adjacency_matrix: np.array,
        seeds: Optional[list[float]] = None,
        names: Optional[list[str]] = None,
        graph_args: Optional[Dict] = None,
        pr_args: Optional[Dict] = None,
    ) -> np.array:
        """
        Base class for algorithm to propagate scores from the adjacency matrix.
        """
        raise NotImplementedError()

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return self.__str__()


class PagerankNX(Ranking):
    """
    Class interfacing NetworkX's Pagerank.
    Initialize with empty constructor if needed.
    """

    @staticmethod
    def propagate(
        adjacency_matrix: np.array,
        seeds: Optional[list[float]] = None,
        names: Optional[list[str]] = None,
        graph_args: Optional[Dict] = None,
        pr_args: Optional[Dict] = None,
    ) -> np.array:
        """
        Propagates scores using NetworkX's Pagerank.
        See NetworkX documentation for input into graph_args and pr_args.
        :param adjacency_matrix: Adjacency matrix
        :type adjacency_matrix: np.array (m x m)
        :param seeds: Metabolite seeds to use as personalization
        :type seeds: List[float]
        :param names: Metabolite names for the metabolites in the graph
        :type names: List[str]
        :param graph_args: Dictionary of arguments passed to networkx DiGraph
        :type graph_args: Optional[dict]
        :param pr_args: Dictionary of arguments passed to networkx Pagerank
        :type pr_args: Optional[dict]
        :return: NumPy array of Pagerank scores
        :rtype: 1-D np.array (m x 1)
        """
        if graph_args is None:
            graph_args = {}
        if pr_args is None:
            pr_args = {}
        graph = nx.DiGraph(adjacency_matrix, **graph_args)
        if isinstance(seeds, list) and len(seeds) > 0:
            pr_args["personalization"] = dict(zip(names, seeds))
            graph = PagerankNX.rename_unnamed_graph(graph, names)
        results = nx.algorithms.link_analysis.pagerank(graph, **pr_args)
        return np.array(list(results.values()))

    @staticmethod
    def rename_unnamed_graph(graph: nx.DiGraph, names: list[str]) -> nx.DiGraph:
        """
        Rename a graph  with unlabeled nodes
        (default node names are integers in range(n_nodes))
        and return a copy labeled with the names given
        :param graph: Graph with unnamed nodes
        :type graph: nx.DiGraph
        :param names: Names to assign to graph nodes
        :type names: list[str]
        :raises ValueError: If number of nodes does not match number of labels
        :return: Identical graph with labeled nodes
        :rtype: nx.DiGraph
        """
        if not len(graph.nodes) == len(names):
            err = "Length of names does not match number of graph nodes"
            logging.error(err)
            raise ValueError(err)
        default_names = range(len(names))
        name_mapping = dict(zip(default_names, names))
        return nx.relabel_nodes(
            graph,
            name_mapping,
        )

    @staticmethod
    def simple_pagerank_nx(adjacency_matrix: np.array) -> np.array:
        """
        Shortcut to run NetworkX's Pagerank (primarily used for testing)
        :param adjacency_matrix: Adjacency matrix
        :type adjacency_matrix: np.array (m x m)
        :return: Pagerank scores
        :rtype: 1-D np.array (m x 1)
        """
        return PagerankNX().propagate(adjacency_matrix)
