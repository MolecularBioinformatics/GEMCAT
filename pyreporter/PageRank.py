#!/usr/bin/python
from sknetwork.ranking import PageRank as PR_SK
import networkx as nx

import numpy as np
import pandas as pd

import pyreporter.utils
from typing import Dict, Tuple, List
import abc


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
    def propagate(A: np.array, graph_args = {}, pr_args = {}) -> np.array:
        """
        Propagates scores using NetworkX's PageRank.
        See NetworkX documentation for input into graph_args and pr_args. 
        :param A: Adjacency matrix
        :type A: np.array (m x m)
        :param graph_args: Dictionary of arguments passed to networkx DiGraph, defaults to {}
        :type graph_args: dict, optional
        :param pr_args: Dictionary of arguments passed to networkx PageRank, defaults to {}
        :type pr_args: dict, optional
        :return: NumPy array of PageRank scores
        :rtype: 1-D np.array (m x 1)
        """
        Ax = nx.DiGraph(A, **graph_args)
        results = nx.algorithms.link_analysis.pagerank(Ax, **pr_args)
        return np.array(list(results.values()))

def run_PR_nx(A: np.array) -> np.array:
    """
    Shortcut to run NetworkX's PageRank (primarily used for testing)
    :param A: Adjacency matrix
    :type A: np.array (m x m)
    :return: n
    :rtype: 1-D np.array (m x 1)
    """
    return PageRankNX().propagate(A)




