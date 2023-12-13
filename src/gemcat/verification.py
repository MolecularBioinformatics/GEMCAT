#!/usr/bin/python

"""
Functions for data verification
"""

import logging
from typing import Union

import pandas as pd


def raise_for_duplicated_index(frame: Union[pd.Series, pd.DataFrame]):
    """
    Check if a DataFrame or Series has duplicated entries in the index.
    :param frame: DataFrame or Series to check
    :type frame: Union[pd.Series, pd.DataFrame]
    :raises ValueError: In case of duplicates
    """
    n_duplicates = frame.index.duplicated(keep="first").sum()
    if n_duplicates > 0:
        err = "Index currently has duplicates"
        logging.error(err)
        raise ValueError(err)


def convert_index_to_str(
    frame: Union[pd.Series, pd.DataFrame]
) -> Union[pd.Series, pd.DataFrame]:
    """
    Checks a DataFrame or Series object's Index dtype.
    If it's not a str, throws a warning and converts it to str.
    :param frame: DataFrame or Series to check
    :type frame: Union[pd.Series, pd.DataFrame]
    :return: DataFrame or Series with Index as a str
    :rtype: Union[pd.Series, pd.DataFrame]
    """
    if not frame.index.dtype == "O":
        logging.warning("Index of DataFrame is not str. Autoconverting to str.")
        frame.index = frame.index.astype(str)
    return frame


def raise_for_non_int_float_64_dtype(series: pd.Series):
    """
    Checks a Series' content for whether it is numeric.
    If it isn't numeric raises a ValueError
    :param series: Series to check
    :type series: pd.Series
    :raises ValueError: ValueError raised if Series content
    isn't NumPy type float64 or int64.
    """
    if not series.dtype == "float64" or series.dtype == "int64":
        err = "Series contains no numerical values"
        logging.error(err)
        raise ValueError(err)
