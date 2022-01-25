import pyreporter.verification as ver
import pandas as pd
from fixtures import *
import pytest

def test_correct_df_index_str():
    testcase = pd.Series({1 : 1})
    expected = pd.Series({'1': 1})
    result = ver.convert_df_index_to_str(testcase)
    assert (result == expected).all()

def test__check_series_content_float():
    testcase = pd.Series({'1': '1'})
    with pytest.raises(ValueError):
        result = ver.raise_for_non_int_float_64_dtype(testcase)

def test_duplicate_indeces_pass():
    testcase = pd.Series({
        'G1': 1.,
        'G2': 2.,
        'G3': 3.,
        })
    ver.raise_for_duplicated_index(testcase)
    assert True

def test_duplicate_indeces_raise():
    testcase = pd.Series({
        'G1': 1.,
        'G2': 2.,
        'G3': 3.,
        })
    testcase = testcase.append(pd.Series({
        'G2': 1.5,
    }))
    with pytest.raises(ValueError):
        ver.raise_for_duplicated_index(testcase)