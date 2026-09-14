#!/usr/bin/python

"""
Collection of smaller utility functions,
mostly related to data processing.
"""

import logging
from typing import List, Tuple, Union

import cobra
import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse as sp


def require_sparse(matrix: object, name: str = "matrix") -> None:
    """
    Raise TypeError unless the argument is a scipy.sparse matrix.

    GEMCAT is sparse-only. Without this check a dense array travels several calls
    deeper and fails on a missing `.data` or `.indptr` instead.
    :param matrix: Object to check
    :type matrix: object
    :param name: Name used in the error message
    :type name: str
    :raises TypeError: If the object is not sparse
    """
    if not sp.issparse(matrix):
        err = (
            f"GEMCAT requires a scipy.sparse matrix for {name}, "
            f"but received {type(matrix).__name__}. "
            f"Convert it with scipy.sparse.csr_array(...)."
        )
        logging.error(err)
        raise TypeError(err)


def prepare_stoich_matrix(stoich_matrix: sp.sparray) -> sp.csc_array:
    """
    Convert a sparse stoichiometric matrix into the form the transforms expect.

    Takes any scipy.sparse format and any numeric dtype. Returns a float64 CSC copy
    with no explicit stored zeros. The input is never modified, so you can transform
    one matrix repeatedly.

    CSC because the next two steps work on columns: scaling by reaction expression,
    and slicing out reversible columns in make_unidirectional.
    :param stoich_matrix: Sparse stoichiometric matrix (m x r)
    :type stoich_matrix: sp.sparray
    :return: CSC copy, float64, no explicit zeros
    :rtype: sp.csc_array
    """
    require_sparse(stoich_matrix, "the stoichiometric matrix")
    # astype copies by default, which protects the caller's matrix.
    working = sp.csc_array(stoich_matrix.astype(np.float64))
    working.eliminate_zeros()
    return working


def as_reaction_vector(expression: npt.ArrayLike, n_reactions: int) -> np.ndarray:
    """
    Flatten an expression vector to a float64 array of length n_reactions.

    Shapes (r,), (1, r) and (r, 1) all mean the same thing: one value per reaction,
    in the column order of the stoichiometric matrix. A wrong length raises instead
    of broadcasting.

    NaN and infinity also raise. Expression scales the stored entries only, so a
    non-finite value would corrupt those and leave the structural zeros at zero: a
    matrix that disagrees with itself. Fill missing genes first; the workflows do.
    :param expression: Expression values shaped (r,), (1, r) or (r, 1)
    :type expression: npt.ArrayLike
    :param n_reactions: Number of reactions the vector must cover
    :type n_reactions: int
    :raises ValueError: On a length mismatch or a non-finite entry
    :return: Flat expression vector
    :rtype: np.ndarray (r,)
    """
    values = np.asarray(expression, dtype=np.float64).ravel()
    if values.size != n_reactions:
        err = (
            f"Expression vector has {values.size} entries "
            f"but the model has {n_reactions} reactions"
        )
        logging.error(err)
        raise ValueError(err)
    if not np.isfinite(values).all():
        err = "Expression vector contains NaN or infinite values"
        logging.error(err)
        raise ValueError(err)
    return values


def scale_columns(matrix: sp.csc_array, factors: np.ndarray) -> sp.csc_array:
    """
    Multiply column j of a CSC matrix by factors[j], in-place destructively.
    This scales reactions by their expression values.
    :param matrix: CSC matrix to scale (mutated in place)
    :type matrix: sp.csc_array
    :param factors: One factor per column
    :type factors: np.ndarray (r,)
    :return: The same matrix, for chaining
    :rtype: sp.csc_array
    """
    matrix.data *= np.repeat(factors, np.diff(matrix.indptr))
    return matrix


def _get_ids(
    iterable: List[Union[cobra.Gene, cobra.Reaction, cobra.Metabolite]]
) -> List[str]:
    """
    Gets a list of IDs from a cobra iterable.
    :param iterable: Iterable for which to get IDs.
    :type iterable: List[Union[cobra.Gene, cobra.Reaction, cobra.Metabolite]]
    :return: List of IDs.
    :rtype: List[str]
    """
    return [g.id for g in iterable]


def _get_n_reactions(stoich_matrix: np.ndarray) -> np.ndarray:
    """
    Returns number of reactions involving each metabolite.
    (number of non-zero entries in a matrix)
    :param stoich_matrix: Stoichiometric matrix (m x r)
    :type stoich_matrix: np.ndarray
    :return: Vector of row-wise sums (total stoichiometries) (m x 1)
    :rtype: np.ndarray
    """
    return np.absolute(np.count_nonzero(stoich_matrix, axis=1))


def _get_total_stoich(stoich_matrix: np.ndarray) -> np.ndarray:
    """
    Returns sum of stoichiometries for each metabolite. (row-wise sums of the matrix).
    Called from within _calc_score_component.
    :param stoich_matrix: Stoichiometric matrix (m x r)
    :type stoich_matrix: np.ndarray
    :return: Vector of row-wise sums (total stoichiometries) (m x 1)
    :rtype: np.ndarray
    """
    return np.absolute(np.sum(stoich_matrix, axis=1))


def split_matrix_pos_neg(matrix: sp.sparray) -> Tuple[sp.sparray, sp.sparray]:
    """
    Split a sparse matrix into its positive and negative parts.
    :param matrix: Sparse matrix to split
    :type matrix: sp.sparray
    :return: Tuple of (positive part, negative part)
    :rtype: Tuple[sp.sparray, sp.sparray]
    """
    require_sparse(matrix)
    epsilon = 0.001

    positive_part = matrix.copy()
    positive_part.data = np.where(matrix.data > epsilon, matrix.data, 0.0)
    positive_part.eliminate_zeros()

    negative_part = matrix.copy()
    negative_part.data = np.where(matrix.data < epsilon, matrix.data, 0.0)
    negative_part.eliminate_zeros()

    return positive_part, negative_part


def annotate_scores(scores: np.ndarray, metabolite_ids: List[str]) -> pd.Series:
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


def get_stoich_matrix_from_cobra(model: cobra.Model) -> sp.csc_array:
    """
    Build the (sparse) stoichiometric matrix of a model.

    Rows follow `model.metabolites`, columns follow `model.reactions`.
    Indices line up with `get_metabolite_ids` and `get_reversibilities` on the same model.
    :param model: Cobra model object for which to get the stoichiometric matrix.
    :type model: cobra.Model
    :return: Stoichiometric matrix (m x r)
    :rtype: sp.csc_array
    """
    metabolite_index = {met.id: i for i, met in enumerate(model.metabolites)}
    rows: List[int] = []
    cols: List[int] = []
    values: List[float] = []
    for column, reaction in enumerate(model.reactions):
        for metabolite, coefficient in reaction.metabolites.items():
            rows.append(metabolite_index[metabolite.id])
            cols.append(column)
            values.append(coefficient)

    matrix = sp.csc_array(
        (values, (rows, cols)),
        shape=(len(metabolite_index), len(model.reactions)),
        dtype=np.float64,
    )
    matrix.eliminate_zeros()
    return matrix


def make_unidirectional(
    stoich_matrix: sp.sparray,
    reversibilities: List[bool],
) -> sp.sparray:
    """
    Split reversible reactions into two opposing one-way reactions.
    Results in an (m x r') matrix with r' = r + the number of reversible reactions.
    :param stoich_matrix: Sparse stoichiometric matrix [m x r]
    :type stoich_matrix: sp.sparray
    :param reversibilities: List of bools whether reactions are reversible
    :type reversibilities: List[bool] [length r]
    :raises TypeError: If a reversibility is not a bool
    :raises ValueError: If the number of reversibilities does not match the columns in the original stoichiometric matrix.
    :return: Unidirectional stoichiometric matrix
    :rtype: sp.sparray [m x r']
    """
    for reversibility_bool in reversibilities:
        if not isinstance(reversibility_bool, bool):
            err = "Bool is expected for reversibility"
            logging.error(err)
            raise TypeError(err)
    require_sparse(stoich_matrix, "the stoichiometric matrix")
    if len(reversibilities) != stoich_matrix.shape[1]:
        err = (
            f"Got {len(reversibilities)} reversibilities "
            f"for {stoich_matrix.shape[1]} reactions"
        )
        logging.error(err)
        raise ValueError(err)

    mask = np.asarray(reversibilities, dtype=bool)
    if not mask.any():
        # Nothing to append
        return stoich_matrix.copy()

    # CSC so the column slice is a contiguous cut of data/indices, not a conversion.
    working = stoich_matrix.tocsc()
    reversible_part = -working[:, np.flatnonzero(mask)]
    return sp.hstack([working, reversible_part], format="csc")


def _get_unidirectional_matrix(model: cobra.Model) -> np.ndarray:
    """
    Takes in a model and returns its stoichiometric matrix
    with reversible reactions separated into
    two different reactions with opposite direction.
    :param model: Model for which to return the stoichiometric matrix.
    :type model: cobra.Model
    :return: Stoichiometric matrix.
    :rtype: np.array (m x r') (where 2r >= r' >= r)
    """
    stoich_matrix = get_stoich_matrix_from_cobra(model)
    reversibilities = [r.reversibility for r in model.reactions]
    return make_unidirectional(stoich_matrix, reversibilities)


def _replace_zeroes(array: np.ndarray) -> np.ndarray:
    """
    Replaces infinity and NaN entries in a matrix with zeroes.
    Called from within _calc_score_component.
    :param array: Array in which to replace values
    :type array: np.ndarray
    :return: Array with entries replaced
    :rtype: np.ndarray
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


def _get_subset_cols(matrix: np.ndarray, indeces: List[int]) -> np.ndarray:
    """
    Get a subset of a matrix by column indeces.
    :param matrix: Matrix to slice.
    :type matrix: np.ndarray
    :param indeces: Indeces of matrix columns to keep.
    :type indeces: List[int]
    :return: Matrix with only given columns included.
    :rtype: np.array (m x |indeces|)
    """
    return matrix[:, indeces]


def _l1_norm(vector: np.ndarray) -> float:
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


def _remove_exchanges(stoich_matrix: np.ndarray, rxn_list: List[str]) -> np.ndarray:
    """
    Remove exchange reactions from a given stoichiometric matrix.
    :param stoich_matrix: Stoichiometric matrix.
    :type stoich_matrix: np.array (m x r)
    :param rxn_list: List of reaction IDs.
    :type rxn_list: List[str]
    :return: Stoichiometric matrix with exchange reactions missing
    :rtype: np.array (m x r' where r' <= r)
    """
    rxn_indeces = _find_indeces(rxn_list)

    return _get_subset_cols(stoich_matrix, rxn_indeces)


def get_reversibilities(model: cobra.Model) -> List[bool]:
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
    return [r.id for r in model.reactions]


def get_metabolite_ids(model: cobra.Model) -> List[str]:
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
    return [m.id for m in model.metabolites]


def make_row_vector(arr: np.ndarray) -> np.ndarray:
    """
    Transform 1D-array into row vector
    :param arr: [description]
    :type arr: np.ndarray
    :raises ValueError: [description]
    :return: [description]
    :rtype: np.ndarray
    """
    return arr.reshape(1, arr.size)


def make_column_vector(arr: np.ndarray) -> np.ndarray:
    """
    Transform 1D-array into column vector
    :param arr: [description]
    :type arr: np.ndarray
    :raises ValueError: [description]
    :return: [description]
    :rtype: np.ndarray
    """
    return arr.reshape(arr.size, 1)


def _is_np_array(arr: object) -> None:
    """
    Throws a TypeError if the object given is not a NumPy array.
    :param arr: Object to check
    :type arr: object
    :raises TypeError: Raised if object type is not np.ndarray
    """
    if not isinstance(arr, np.ndarray):
        received = type(arr)
        err = f"Expected a NumPy array but received {received}"
        logging.error(err)
        raise TypeError(err)


def _check_array_shape(arr: np.ndarray, target: np.ndarray) -> None:
    """
    Raise ValueError unless two arrays have the same shape.
    :param arr: Array to check
    :type arr: np.ndarray
    :param target: Array whose shape arr must match
    :type target: np.ndarray
    :raises ValueError: Raised if shape of the two arrays doesn't match.
    """
    if not arr.shape == target.shape:
        msg = f"Array shape needs to be {target.shape} but is {arr.shape}"
        logging.error(msg)
        raise ValueError(msg)


def _is_all_ones(arr: np.ndarray) -> bool:
    """
    Returns true if an array is all ones.
    :param arr: Array to check
    :type arr: np.ndarray
    :return: True if array is all ones
    :rtype: bool
    """
    ones = np.ones(arr.shape)
    return np.allclose(arr, ones)


def geometric_mean(*numbers: float) -> float:
    """
    Calculates the geometric mean for a number of ints or floats.
    :param numbers: Variable number of floats or ints (vararg)
    :type numbers: float
    :raises ValueError: If called with no numbers
    :return: Geometric mean
    :rtype: float
    """
    n_numbers = len(numbers)
    if n_numbers == 0:
        err = "Cannot calculate the geometric mean of an empty set of numbers"
        logging.error(err)
        raise ValueError(err)
    as_floats = [float(i) for i in numbers]
    prod = multiply(as_floats)
    return float(prod ** (1 / n_numbers))


def multiply(*numbers: npt.ArrayLike) -> np.number:
    """
    Convenience function to get the product of a list of items.
    Takes loose scalars or one sequence; geometric_mean passes a list.
    :param numbers: Variable number of floats, ints, or one sequence of them (vararg)
    :type numbers: npt.ArrayLike
    :return: Product of the numbers
    :rtype: np.number
    """
    return np.prod(np.asarray(numbers, dtype=np.float64))
