#!/usr/bin/python

"""
Algorithms relating to calculation of the adjacency matrix
"""

import abc

import numpy as np

from . import utils


class AdjacencyTransformation(abc.ABC):
    """
    Abstract base class defining adjacency transformation interface do NOT use.
    """

    @staticmethod
    @abc.abstractmethod
    def transform(
        stoich_matrix: np.array,
        reversibilities: list[bool],
        expression: np.array,
    ) -> np.array:
        """
        Calculates the (Pagerank) adjacency matrix for a given stoichiometric matrix.
        :param stoich_matrix: Stoichiometric matrix.
        :type stoich_matrix: np.array (m x r)
        :param reversibilities: Reaction reversibilities
        :type reversibilities: list[bool]
        :param expression: Reaction expression vector
        :type expression: np.array
        :return: Pagerank adjacency matrix.
        :rtype: np.array (m x m)
        """
        raise NotImplementedError()


class ATFullStoich(AdjacencyTransformation):
    """
    Transform S into A with normalization based on rows, stoichiometry in both
    products and educts considered.
    """

    @staticmethod
    def transform(
        stoich_matrix: np.array,
        reversibilities: list[bool],
        expression: np.array,
    ) -> np.array:
        """
        Calculates the (Pagerank) adjacency matrix for a given stoichiometric matrix.
        :param stoich_matrix: Stoichiometric matrix.
        :type stoich_matrix: np.array (m x r)
        :param reversibilities: Reaction reversibilities
        :type reversibilities: list[bool]
        :param expression: Reaction expression vector
        :type expression: np.array
        :return: Pagerank adjacency matrix.
        :rtype: np.array (m x m)
        """
        stoich_matrix = stoich_matrix * expression
        stoich_matrix = utils.make_unidirectional(stoich_matrix, reversibilities)
        epsilon = 10 ** (-20)
        positive_part, negative_part = utils.split_matrix_pos_neg(stoich_matrix)

        adjacencies = np.abs(negative_part) @ np.transpose(positive_part)

        col_sum = utils.make_column_vector(adjacencies.sum(axis=1) + epsilon)
        adjacencies = np.divide(adjacencies, col_sum)

        return adjacencies


class ATHalfStoich(AdjacencyTransformation):
    """
    Transform S into A with normalization based on rows,
    stoichiometry in products only considered.
    """

    @staticmethod
    def transform(
        stoich_matrix: np.array,
        reversibilities: list[bool],
        expression: np.array,
    ) -> np.array:
        """
        Calculates the (Pagerank) adjacency matrix for a given stoichiometric matrix.
        :param stoich_matrix: Stoichiometric matrix.
        :type stoich_matrix: np.array (m x r)
        :param reversibilities: Reaction reversibilities
        :type reversibilities: list[bool]
        :param expression: Expression data vector (1 x r)
        :type expression: np.array
        :return: Pagerank adjacency matrix.
        :rtype: np.array (m x m)
        """
        stoich_matrix = stoich_matrix * expression
        stoich_matrix = utils.make_unidirectional(stoich_matrix, reversibilities)
        epsilon = 10 ** (-20)
        positive_part, negative_part = utils.split_matrix_pos_neg(stoich_matrix)

        outgoing_ones = np.divide(negative_part, (negative_part + epsilon))
        adjacencies = outgoing_ones @ np.transpose(positive_part)

        col_sum = utils.make_column_vector(adjacencies.sum(axis=1) + epsilon)
        adjacencies = np.divide(adjacencies, col_sum)

        return adjacencies


class ATPureAdjacency(AdjacencyTransformation):
    """
    Transform S into A using pure adjacency not regarding stoichiometry.
    """

    @staticmethod
    def transform(
        stoich_matrix: np.array,
        reversibilities: list[bool],
        expression: np.array,
    ) -> np.array:
        """
        Calculates the (Pagerank) adjacency matrix for a given stoichiometric matrix.
        :param stoich_matrix: Stoichiometric matrix.
        :type stoich_matrix: np.array (m x r)
        :param reversibilities: Reaction reversibilities
        :type reversibilities: list[bool]
        :param expression: Expression data array (1 x r)
        :type expression: np.array
        :return: Pagerank adjacency matrix.
        :rtype: np.array (m x m)
        """
        epsilon = 10 ** (-20)
        stoich_matrix = stoich_matrix / np.abs(stoich_matrix + epsilon)
        stoich_matrix = stoich_matrix * utils.make_row_vector(expression)
        stoich_matrix = utils.make_unidirectional(stoich_matrix, reversibilities)
        positive_part, negative_part = utils.split_matrix_pos_neg(stoich_matrix)

        outgoing_ones = np.divide(negative_part, (negative_part + epsilon))
        adjacencies = outgoing_ones @ np.transpose(positive_part)

        col_sum = utils.make_column_vector(adjacencies.sum(axis=1) + epsilon)
        adjacencies = np.divide(adjacencies, col_sum)

        return adjacencies


def run_adjacencies_normalize(
    stoich_matrix: np.array,
    reversibilities: list[bool],
    expression_vector: np.array,
    adjacency_transformation: AdjacencyTransformation,
) -> np.array:
    """
    Shortcut function for testing calculation of A with row-based normalization.
    :param stoich_matrix: stoichiometric matrix
    :type stoich_matrix: np.array (m x r)
    :param reversibilities: list of reaction reversibilities
    :type reversibilities: list[bool]
    :param expression_vector: Expression data array
    :type expression_vector: np.array (1 x r)
    :param adjacency_transformation: Adjacency transformation algorithm
    :type adjacency_transformation: AdjacencyTransformation
    :return: Adjacency matrix A
    :rtype: np.array (m x m)
    """
    adjacency = adjacency_transformation()
    return adjacency.transform(stoich_matrix, reversibilities, expression_vector)
