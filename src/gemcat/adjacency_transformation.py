#!/usr/bin/python

"""
Algorithms relating to calculation of the adjacency matrix
"""

import abc
from typing import Type

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from . import utils

# Guards division by a zero row sum, for a metabolite with no outgoing edges.
EPSILON = 10 ** (-20)


class AdjacencyTransformation(abc.ABC):
    """
    Turns a stoichiometric matrix into a metabolite adjacency matrix.
    Abstract base class to define the interface.

    Every current implementation runs the same five steps:

      1. weight each reaction column by its expression value
      2. split reversible reactions into two opposing one-way reactions
      3. separate each reaction into educts and products
      4. link each educt to each product
      5. scale each metabolite's outgoing edges to sum to 1

    The result A is (m x m). A[i, j] is the weighted edge from metabolite i to metabolite j.
    Rows sum to 1.0, or to 0.0 for a metabolite with no outgoing edges.

    Current implementations differ by how much stoichiometry shapes the ranking.
    ATPureAdjacency ignores coefficients and is the default; ATHalfStoich keeps them
    for products only; ATFullStoich keeps them on educts and products.
    """

    @staticmethod
    @abc.abstractmethod
    def transform(
        stoich_matrix: sp.sparray,
        reversibilities: list[bool],
        expression: npt.ArrayLike,
    ) -> sp.csr_array:
        """
        Calculates the (Pagerank) adjacency matrix for a given stoichiometric matrix.
        :param stoich_matrix: Sparse stoichiometric matrix (m x r)
        :type stoich_matrix: sp.sparray
        :param reversibilities: Reaction reversibilities
        :type reversibilities: list[bool]
        :param expression: Reaction expression, shaped (r,), (1, r) or (r, 1)
        :type expression: npt.ArrayLike
        :return: Adjacency matrix; A[i, j] is the edge weight from metabolite i to j
        :rtype: sp.csr_array (m x m)
        """
        raise NotImplementedError()


def _finalize_adjacencies(
    outgoing: sp.sparray, positive_part: sp.sparray
) -> sp.csr_array:
    """
    Shared post-processing of all transforms: pair educts with products, normalize rows.
    Dividing each row by its total turns a metabolite's outgoing weights into a probability
    distribution for PageRank.

    Rows sum to 1.0 (0.0 for a metabolite with no outgoing edges).
    :param outgoing: Left matmul operand (m x r')
    :type outgoing: sp.sparray
    :param positive_part: Right operand, transposed internally (m x r')
    :type positive_part: sp.sparray
    :raises TypeError: If the row sums are not 1-D, i.e. a legacy spmatrix leaked in
    :return: Row-normalized adjacency matrix
    :rtype: sp.csr_array
    """
    # positive_part.T is a zero-copy view over the same buffers. Do not mutate.
    adjacencies = sp.csr_array(outgoing.tocsr() @ positive_part.T)
    # remove explicit zeros to avoid zero-weight edges in the graph
    adjacencies.eliminate_zeros()

    row_sums = np.asarray(adjacencies.sum(axis=1)).ravel()
    if row_sums.shape != (adjacencies.shape[0],):
        # Only a legacy spmatrix gets here: sum(axis=1) returns a 2-D np.matrix,
        # which would broadcast into a dense (m x m) array. Fail instead.
        raise TypeError(
            f"Expected 1-D row sums of length {adjacencies.shape[0]}, "
            f"got {row_sums.shape}"
        )

    # np.repeat expands the row sums into a per-entry divisor.
    # Divide rather than multiply by a reciprocal to avoid rounding twice.
    adjacencies.data /= np.repeat(row_sums + EPSILON, np.diff(adjacencies.indptr))
    return adjacencies


class ATFullStoich(AdjacencyTransformation):
    """
    Adjacency weighted by the stoichiometry of both educts and products.

    The edge from educt i to product j carries |coefficient of i| * coefficient of
    j. A reaction consuming 2 ATP to make 3 NADH therefore weighs six times a 1:1
    conversion.

    See AdjacencyTransformation for the shared steps.
    """

    @staticmethod
    def transform(
        stoich_matrix: sp.sparray,
        reversibilities: list[bool],
        expression: npt.ArrayLike,
    ) -> sp.csr_array:
        """
        Calculates the (Pagerank) adjacency matrix for a given stoichiometric matrix.
        :param stoich_matrix: Sparse stoichiometric matrix (m x r)
        :type stoich_matrix: sp.sparray
        :param reversibilities: Reaction reversibilities
        :type reversibilities: list[bool]
        :param expression: Reaction expression, shaped (r,), (1, r) or (r, 1)
        :type expression: npt.ArrayLike
        :return: Adjacency matrix; A[i, j] is the edge weight from metabolite i to j
        :rtype: sp.csr_array (m x m)
        """
        working = utils.prepare_stoich_matrix(stoich_matrix)
        expression = utils.as_reaction_vector(expression, working.shape[1])

        utils.scale_columns(working, expression)
        working = utils.make_unidirectional(working, reversibilities)
        positive_part, negative_part = utils.split_matrix_pos_neg(working)

        negative_part.data = np.abs(negative_part.data)

        return _finalize_adjacencies(negative_part, positive_part)


class ATHalfStoich(AdjacencyTransformation):
    """
    Adjacency weighted by product stoichiometry only.

    Every educt counts the same whatever its coefficient; products keep theirs.

    See AdjacencyTransformation for the shared steps.
    """

    @staticmethod
    def transform(
        stoich_matrix: sp.sparray,
        reversibilities: list[bool],
        expression: npt.ArrayLike,
    ) -> sp.csr_array:
        """
        Calculates the (Pagerank) adjacency matrix for a given stoichiometric matrix.
        :param stoich_matrix: Sparse stoichiometric matrix (m x r)
        :type stoich_matrix: sp.sparray
        :param reversibilities: Reaction reversibilities
        :type reversibilities: list[bool]
        :param expression: Reaction expression, shaped (r,), (1, r) or (r, 1)
        :type expression: npt.ArrayLike
        :return: Adjacency matrix; A[i, j] is the edge weight from metabolite i to j
        :rtype: sp.csr_array (m x m)
        """
        working = utils.prepare_stoich_matrix(stoich_matrix)
        expression = utils.as_reaction_vector(expression, working.shape[1])

        utils.scale_columns(working, expression)
        working = utils.make_unidirectional(working, reversibilities)
        positive_part, negative_part = utils.split_matrix_pos_neg(working)

        # Make stoichiometry 1 or 0 everywhere
        negative_part.data /= negative_part.data + EPSILON

        return _finalize_adjacencies(negative_part, positive_part)


class ATPureAdjacency(AdjacencyTransformation):
    """
    Adjacency that ignores stoichiometry. GEMCAT's default.

    Stoichiometric coefficients are stripped; expression and topology set every
    edge weight.

    See AdjacencyTransformation for the shared steps.
    """

    @staticmethod
    def transform(
        stoich_matrix: sp.sparray,
        reversibilities: list[bool],
        expression: npt.ArrayLike,
    ) -> sp.csr_array:
        """
        Calculates the (Pagerank) adjacency matrix for a given stoichiometric matrix.
        :param stoich_matrix: Sparse stoichiometric matrix (m x r)
        :type stoich_matrix: sp.sparray
        :param reversibilities: Reaction reversibilities
        :type reversibilities: list[bool]
        :param expression: Reaction expression, shaped (r,), (1, r) or (r, 1)
        :type expression: npt.ArrayLike
        :return: Adjacency matrix; A[i, j] is the edge weight from metabolite i to j
        :rtype: sp.csr_array (m x m)
        """
        working = utils.prepare_stoich_matrix(stoich_matrix)
        expression = utils.as_reaction_vector(expression, working.shape[1])

        # Make stoichiometry one, zero, or minus one everywhere
        working.data /= np.abs(working.data + EPSILON)
        utils.scale_columns(working, expression)
        working = utils.make_unidirectional(working, reversibilities)
        positive_part, negative_part = utils.split_matrix_pos_neg(working)

        # make stoichiometry one or zero everywhere in educts
        negative_part.data /= negative_part.data + EPSILON

        return _finalize_adjacencies(negative_part, positive_part)


def run_adjacencies_normalize(
    stoich_matrix: sp.sparray,
    reversibilities: list[bool],
    expression_vector: npt.ArrayLike,
    adjacency_transformation: Type[AdjacencyTransformation],
) -> sp.csr_array:
    """
    Shortcut function for testing calculation of A with row-based normalization.
    :param stoich_matrix: Sparse stoichiometric matrix (m x r)
    :type stoich_matrix: sp.sparray
    :param reversibilities: list of reaction reversibilities
    :type reversibilities: list[bool]
    :param expression_vector: Reaction expression, shaped (r,), (1, r) or (r, 1)
    :type expression_vector: npt.ArrayLike
    :param adjacency_transformation: Adjacency transformation class, not an
        instance; this function instantiates it
    :type adjacency_transformation: Type[AdjacencyTransformation]
    :return: Adjacency matrix A
    :rtype: sp.csr_array (m x m)
    """
    adjacency = adjacency_transformation()
    return adjacency.transform(stoich_matrix, reversibilities, expression_vector)
