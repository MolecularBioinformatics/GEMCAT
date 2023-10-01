#!/usr/bin/python
import abc
import logging
from typing import Dict, Optional

import networkx as nx
import numpy as np


class Ranking(abc.ABC):
    """
    Abstract base class for Ranking methods. Do NOT use.
    """

    @abc.abstractmethod
    def propagate(self, graph_args: Dict, pr_args: Dict) -> np.array:
        pass


class PageRankNX(Ranking):
    """
    Class interfacing NetworkX's PageRank.
    Initialize with empty constructor if needed.
    """

    @staticmethod
    def propagate(
        A: np.array,
        seeds: Optional[list[float]] = None,
        names: Optional[list[str]] = None,
        graph_args: Optional[Dict] = None,
        pr_args: Optional[Dict] = None,
    ) -> np.array:
        """
        Propagates scores using NetworkX's PageRank.
        See NetworkX documentation for input into graph_args and pr_args.
        :param A: Adjacency matrix
        :type A: np.array (m x m)
        :param seeds: Metabolite seeds to use as personalization
        :type seeds: List[float]
        :param names: Metabolite names for the metabolites in the graph
        :type names: List[str]
        :param graph_args: Dictionary of arguments passed to networkx DiGraph, defaults to None
        :type graph_args: Optional[dict]
        :param pr_args: Dictionary of arguments passed to networkx PageRank, defaults to None
        :type pr_args: Optional[dict]
        :return: NumPy array of PageRank scores
        :rtype: 1-D np.array (m x 1)
        """
        if graph_args is None:
            graph_args = {}
        if pr_args is None:
            pr_args = {}
        Ax = nx.DiGraph(A, **graph_args)
        if isinstance(seeds, list) and len(seeds) > 0:
            pr_args["personalization"] = dict(zip(names, seeds))
            Ax = rename_unnamed_graph(Ax, names)
        results = nx.algorithms.link_analysis.pagerank(Ax, **pr_args)
        return np.array(list(results.values()))


def rename_unnamed_graph(G: nx.DiGraph, names: list[str]) -> nx.DiGraph:
    """_summary_
    Rename a graph  with unlabeled nodes
    (default node names are integers in range(n_nodes))
    and return a copy labeled with the names given
    :param G: Graph with unnamed nodes
    :type G: nx.DiGraph
    :param names: Names to assign to graph nodes
    :type names: list[str]
    :raises ValueError: If number of nodes does not match number of labels
    :return: Identical graph with labeled nodes
    :rtype: nx.DiGraph
    """
    if not len(G.nodes) == len(names):
        err = "Length of names does not match number of graph nodes"
        logging.error(err)
        raise ValueError(err)
    default_names = range(len(names))
    name_mapping = dict(zip(default_names, names))
    return nx.relabel_nodes(
        G,
        name_mapping,
    )


def run_PR_nx(A: np.array) -> np.array:
    """
    Shortcut to run NetworkX's PageRank (primarily used for testing)
    :param A: Adjacency matrix
    :type A: np.array (m x m)
    :return: n
    :rtype: 1-D np.array (m x 1)
    """
    return PageRankNX().propagate(A)
