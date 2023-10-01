import abc
from typing import List

import numpy as np

from . import utils


class AdjacencyTransformation(abc.ABC):
    """
    Abstract base class defining adjacency transformation interface do NOT use.
    """

    @abc.abstractmethod
    def transform(S: np.array) -> np.array:
        pass


class ATFullStoich(AdjacencyTransformation):
    """
    Transform S into A with normalization based on rows, stoichiometry in both
    products and educts considered.
    """

    @staticmethod
    def transform(
        S: np.array,
        reversibilities: List[bool],
        expression: np.array,
    ) -> np.array:
        """
        Calculates the (PageRank) adjacency matrix for a given stoichiometric matrix.
        :param S: Stoichiometric matrix.
        :type S: np.array (m x r)
        :return: PageRank adjacency matrix.
        :rtype: np.array (m x m)
        """
        S = S * expression
        S = utils._make_unidirectional(S, reversibilities)
        c = 10 ** (-20)
        S_plus, S_minus = utils._split_matrix_pos_neg(S)

        A = np.abs(S_minus) @ np.transpose(S_plus)

        n = utils._make_column_vector(A.sum(axis=1) + c)
        A = np.divide(A, n)

        return A


class ATHalfStoich(AdjacencyTransformation):
    """
    Transform S into A with normalization based on rows,
    stoichiometry in products only considered.
    """

    @staticmethod
    def transform(
        S: np.array,
        reversibilities: List[bool],
        expression: np.array,
    ) -> np.array:
        """
        Calculates the (PageRank) adjacency matrix for a given stoichiometric matrix.
        :param S: Stoichiometric matrix.
        :type S: np.array (m x r)
        :return: PageRank adjacency matrix.
        :rtype: np.array (m x m)
        """
        S = S * expression
        S = utils._make_unidirectional(S, reversibilities)
        c = 10 ** (-20)
        S_plus, S_minus = utils._split_matrix_pos_neg(S)

        S_out = np.divide(S_minus, (S_minus + c))
        A = S_out @ np.transpose(S_plus)

        n = utils._make_column_vector(A.sum(axis=1) + c)
        A = np.divide(A, n)

        return A


class ATPureAdjacency(AdjacencyTransformation):
    """
    Transform S into A using pure adjacency not regarding stoichiometry.
    """

    @staticmethod
    def transform(
        S: np.array,
        reversibilities: List[bool],
        expression: np.array,
    ) -> np.array:
        """
        Calculates the (PageRank) adjacency matrix for a given stoichiometric matrix.
        :param S: Stoichiometric matrix.
        :type S: np.array (m x r)
        :return: PageRank adjacency matrix.
        :rtype: np.array (m x m)
        """
        c = 10 ** (-20)
        S = S / np.abs(S + c)
        S = S * utils._make_row_vector(expression)
        S = utils._make_unidirectional(S, reversibilities)
        S_plus, S_minus = utils._split_matrix_pos_neg(S)

        S_out = np.divide(S_minus, (S_minus + c))
        A = S_out @ np.transpose(S_plus)

        n = utils._make_column_vector(A.sum(axis=1) + c)
        A = np.divide(A, n)

        return A


def run_AT_normalize(
    S: np.array,
    reversibilities: List[bool],
    expression_vector: np.array,
    AT: AdjacencyTransformation,
) -> np.array:
    """
    Shortcut function for testing calculation of A with row-based normalization.
    :param S: stoichiometric matrix
    :type S: np.array (m x r)
    :return: Adjacency matrix A
    :rtype: np.array (m x m)
    """
    at = AT()
    return at.transform(S, reversibilities, expression_vector)
