import logging
from typing import Union

import pandas as pd


def raise_for_duplicated_index(df: Union[pd.Series, pd.DataFrame]):
    n_duplicates = df.index.duplicated(keep="first").sum()
    if n_duplicates > 0:
        err = "Index currently has duplicates"
        logging.error(err)
        raise ValueError(err)


def convert_df_index_to_str(
    df: Union[pd.Series, pd.DataFrame]
) -> Union[pd.Series, pd.DataFrame]:
    """
    Checks a DataFrame or Series object's Index dtype. If it's not a str, throws a warning and converts it to str.
    :param df: DataFrame or Series to check
    :type df: Union[pd.Series, pd.DataFrame]
    :return: DataFrame or Series with Index as a str
    :rtype: Union[pd.Series, pd.DataFrame]
    """
    if not df.index.dtype == "O":
        logging.warning("Index of DataFrame is not str. Autoconverting to str.")
        df.index = df.index.astype(str)
    return df


def raise_for_non_int_float_64_dtype(ser: pd.Series):
    """
    Checks a Series' content for whether it is numeric. If it isn't numeric raises a ValueError
    :param ser: Series to check
    :type ser: pd.Series
    :raises ValueError: ValueError raised if Series content isn't NumPy type float64 or int64.
    """
    if not ser.dtype == "float64" or ser.dtype == "int64":
        err = "Series contains no numerical values"
        logging.error(err)
        raise ValueError(err)
