#!/usr/bin/python

import logging
import math
from functools import reduce
from typing import Iterable, List, Tuple, Union

import cobra
import numpy as np
import pandas as pd


def _get_ids(
    iterable: List[Union[cobra.Gene, cobra.Reaction, cobra.Metabolite]]
) -> List[str]:
    """
    Gets a list of IDs from a cobra iterable.
    :param genes: Iterable for which to get IDs.
    :type genes: List[Union[cobra.Gene, cobra.Reaction, cobra.Metabolite]]
    :return: List of IDs.
    :rtype: List[str]
    """
    return [g.id for g in iterable]


def _get_n_reactions(s: np.array) -> np.array:
    """
    Returns number of reactions involving each metabolite. (number of non-zero entries in a matrix)
    :param s: Stoichiometric matrix (m x r)
    :type s: np.array
    :return: Vector of row-wise sums (total stoichiometries) (m x 1)
    :rtype: np.array
    """
    return np.absolute(np.count_nonzero(s, axis=1))


def _get_total_stoich(s: np.array) -> np.array:
    """
    Returns sum of stoichiometries for each metabolite. (row-wise sums of the matrix).
    Called from within _calc_score_component.
    :param s: Stoichiometric matrix (m x r)
    :type s: np.array
    :return: Vector of row-wise sums (total stoichiometries) (m x 1)
    :rtype: np.array
    """
    return np.absolute(np.sum(s, axis=1))


def _split_matrix_pos_neg(s: np.array) -> Tuple[np.array, np.array]:
    """
    Splits an array into two arrays, one containing all positive entries, one containing all negative entries.
    Called from within _calc_met_score.
    :param s: Array to split
    :type s: np.array
    :return: Tuple of (postive array, negative array)
    :rtype: Tuple[np.array, np.array]
    """
    s_plus = s * (s > 0.001)
    s_minus = s * (s < 0.001)
    return s_plus, s_minus


def _annotate(scores: np.array, metabolite_ids: List[str]) -> pd.Series:
    """
    Create a pandas Series matching metabolite scores with their IDs.
    :param scores: Metabolite scores.
    :type scores: np.array (m x 1)
    :param metabolite_ids: List of metabolite IDs.
    :type metabolite_ids: List[str]
    :return: Series matching scores and IDs.
    :rtype: pd.Series (m x 1)
    """
    return pd.Series(scores, index=metabolite_ids)


def _get_stoich_matrix(model: cobra.Model) -> np.array:
    """
    Returns the stoichiometric matrix of a model.
    Called from within get_initial_scores.
    :param model: Cobra model object for which to get the stoichiometric matrix.
    :type model: cobra.Model
    :return: Stoichiometric matrix (m x r)
    :rtype: np.array
    """
    return cobra.util.array.create_stoichiometric_matrix(
        model, array_type="dense", dtype=float
    )


def _make_unidirectional(
    S: np.array,
    reversibilities: List[bool],
) -> np.array:
    """
    Takes in a stoichiometric matrix and a list of reversibilities,
    then adds the reverse of all the reversible reactions to the matrix.
    :param S: Stoichiometric matrix
    :type S: np.array [m x r]
    :param reversibilities: List of bools whether reactions are reversible
    :type reversibilities: List[bool] [length r]
    :return: Unidirectional stoichiometric matrix
    :rtype: np.array [m x r']
    """
    for r in reversibilities:
        if not isinstance(r, bool):
            err = "Bool is expected for reversibility"
            logging.error(err)
            raise TypeError(err)
    S_rev = S[:, reversibilities]
    S_rev = -1.0 * S_rev
    return np.append(S, S_rev, axis=1)


def _get_unidirectional_matrix(model: cobra.Model) -> np.array:
    """
    Takes in a model and returns its stoichiometric matrix
    with reversible reactions separated into
    two different reactions with opposite direction.
    :param model: Model for which to return the stoichiometric matrix.
    :type model: cobra.Model
    :return: Stoichiometric matrix.
    :rtype: np.array (m x r') (where 2r >= r' >= r)
    """
    S = _get_stoich_matrix(model)
    reversibilities = [r.reversibility for r in model.reactions]
    return _make_unidirectional(S, reversibilities)


def _replace_zeroes(array: np.array) -> np.array:
    """
    Replaces infinity and NaN entries in a matrix with zeroes.
    Called from within _calc_score_component.
    :param array: Array in which to replace values
    :type array: np.array
    :return: Array with entries replaced
    :rtype: np.array
    """
    array[array == np.inf] = 0.0
    array[array == -np.inf] = 0.0
    array[np.isnan(array)] = 0.0
    return array


def _calc_zscore(series: pd.Series) -> pd.Series:
    """
    Calculate z-Score of a Pandas Series.
    :param series: Series to calculate z-Score of.
    :type series: pd.Series
    :return: Pandas Series of z-Score.
    :rtype: pd.Series
    """
    return (series - series.mean()) / series.std()


def _scale(series: pd.Series) -> pd.Series:
    """
    Scale a Pandas Series .
    :param series: Series of scores to scale.
    :type series: pd.Series
    :return: Pandas Series of scaled scores.
    :rtype: pd.Series
    """
    # return series / max(abs(series.min()), series.max())
    return series / series.sum()


def _find_indeces(rxn_list: List[str]) -> List[int]:
    """
    Out of a list of reaction strings, find the indeces of the non-exchange reactions.
    :param rxn_list: List of reaction indeces
    :type rxn_list: List[str]
    :return: List of indeces of non-exchange reactions
    :rtype: List[int]
    """
    enum = enumerate(rxn_list)
    filtered = [count for (count, tag) in enum if not _is_exchange(tag)]

    return filtered


def _is_exchange(tag: str) -> bool:
    """
    Use the reaction ID to determine whether it is an exchange reaction.
    :param tag: Reaction ID
    :type tag: str
    :return: True/False whether reaction is an exchange reaction.
    :rtype: bool
    """
    exchange_prefixes = ["OF_", "EX_"]
    for prefix in exchange_prefixes:
        if tag.startswith(prefix):
            return True

    return False


def _get_subset_cols(s: np.array, indeces: List[int]) -> np.array:
    """
    Get a subset of a matrix by column indeces.
    :param s: Matrix to slice.
    :type s: np.array
    :param indeces: Indeces of matrix columns to keep.
    :type indeces: List[int]
    :return: Matrix with only given columns included.
    :rtype: np.array (m x |indeces|)
    """

    return s[:, indeces]


def _l1_norm(vector: np.array) -> float:
    """
    Returns the L1-Norm (Manhattan distance) of a NumPy array.
    :param vector: Vector of which to calculate the L1-Norm
    :type vector: np.array (m x 1)
    :return: L1-norm of vector
    :rtype: float
    """
    if vector.size == 0:
        err = "Cannot calculate the l1-norm of an empty vector"
        logging.error(err)
        raise ValueError(err)
    return np.sum(np.abs(vector))


def _remove_exchanges(s: np.array, rxn_list: List[str]) -> np.array:
    """
    Remove exchange reactions from a given stoichiometric matrix.
    :param s: Stoichiometric matrix.
    :type s: np.array (m x r)
    :param rxn_list: List of reaction IDs.
    :type rxn_list: List[str]
    :return: Stoichiometric matrix with exchange reactions missing
    :rtype: np.array (m x r' where r' <= r)
    """
    rxn_indeces = _find_indeces(rxn_list)

    return _get_subset_cols(s, rxn_indeces)


def _get_reversibilities(model: cobra.Model) -> List[bool]:
    """
    Return a list of reversibilities for the model.
    :param model: Model from which to extract reversibilities
    :type model: cobra.Model
    :return: List of reversibilities. True for reversible reactions.
    :rtype: List[bool]
    """
    return [r.reversibility for r in model.reactions]


def _get_reaction_ids(model: cobra.Model) -> List[str]:
    """
    Returns the list of reaction IDs from a given model.
    :param model: Model object
    :type model: cobra.Model
    :return: List of reaction IDs in the model
    :rtype: List[str]
    """
    if not isinstance(model, cobra.Model):
        err = "CobraPy model required to extract reaction IDs"
        logging.error(err)
        raise TypeError(err)
    if len(model.reactions) == 0:
        err = "The COBRA model contains no reactions"
        logging.error(err)
        raise ValueError(err)
    r_ids = [r.id for r in model.reactions]
    return r_ids


def _get_metabolite_ids(model: cobra.Model) -> List[str]:
    """
    Returns the list of metabolite IDs from a given model.
    :param model: Model object
    :type model: cobra.Model
    :return: List of metabolite IDs in the model
    :rtype: List[str]
    """
    if not isinstance(model, cobra.Model):
        err = "CobraPy model required to extract metabolite IDs"
        logging.error(err)
        raise TypeError(err)
    if len(model.metabolites) == 0:
        err = "The COBRA model contains no metabolites"
        logging.error(err)
        raise ValueError(err)
    m_ids = [m.id for m in model.metabolites]
    return m_ids


def _make_row_vector(arr: np.array) -> np.array:
    """
    Transform 1D-array into row vector
    :param arr: [description]
    :type arr: np.array
    :raises ValueError: [description]
    :return: [description]
    :rtype: np.array
    """
    return arr.reshape(1, arr.size)


def _make_column_vector(arr: np.array) -> np.array:
    """
    Transform 1D-array into column vector
    :param arr: [description]
    :type arr: np.array
    :raises ValueError: [description]
    :return: [description]
    :rtype: np.array
    """
    return arr.reshape(arr.size, 1)


def _is_np_array(arr):
    """
    Throws a TypeError if the object given is not a NumPy array.
    :param arr: Object to check
    :type arr: Object
    :raises TypeError: Raised if object type is not np.array
    """
    if not isinstance(arr, np.ndarray):
        received = type(arr)
        err = f"Expected a NumPy array but received {received}"
        logging.error(err)
        raise TypeError(err)


def _check_array_shape(arr, target):
    """
    :param arr: Array to check
    type arr: np.array
    :raises ValueError: Raised if shape of the two arrays doesn't match.
    """
    if not arr.shape == target.shape:
        msg = f"Array shape needs to be {target.shape} but is {arr.shape}"
        logging.error(msg)
        raise ValueError(msg)


def _is_all_ones(arr: np.array) -> bool:
    """
    Returns true if an array is all ones.
    :param arr: Array to check
    :type arr: np.array
    :return: True if array is all ones
    :rtype: bool
    """
    ones = np.ones(arr.shape)
    return np.allclose(arr, ones)


def geometric_mean(*numbers):
    """
    Calculates the geometric mean for a number of ints or floats.
    :return: Geometric mean
    :rtype: float
    """
    n = len(numbers)
    if n == 0:
        err = "Cannot calculate the geometric mean of an empty set of numbers"
        logging.error(err)
        raise ValueError(err)
    numbers = [float(i) for i in numbers]
    prod = multiply(numbers)
    return prod ** (1 / n)


def multiply(*numbers):
    return np.prod(numbers)
