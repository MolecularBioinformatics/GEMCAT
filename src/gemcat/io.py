"""
Input and output functionalities of the framework
"""

import logging
import pickle
from pathlib import Path
from typing import Optional, Union

import cobra
import pandas as pd

from . import utils
from .model import Model


def convert_cobra_model(cobra_model: cobra.Model) -> Model:
    """
    Converts a COBRA model into a Pagerank-based Model.
    :param cobra_model: COBRA model object
    :type cobra_model: cobra.Model
    :return: Page-Rank based Model
    :rtype: Model
    """
    if not isinstance(cobra_model, cobra.Model):
        received = type(cobra_model)
        err = f"Expected a COBRA Model but received {received}"
        logging.error(err)
        raise TypeError(err)
    stoich_matrix = utils.get_stoich_matrix_from_cobra(cobra_model)
    rev = utils.get_reversibilities(cobra_model)
    metabolites = utils.get_metabolite_ids(cobra_model)
    return Model(stoich_matrix, metabolites, rev)


def load_json_cobra(json_file: Union[str, Path]) -> cobra.Model:
    """
    Loads a Pagerank-based Model from a COBRA json file.
    :param json_file: Path to json model file
    :type json_file: Union[str, Path]
    :return: Model object created from the COBRA model
    :rtype: Model
    """
    cobra_model = cobra.io.load_json_model(str(json_file))
    return cobra_model


def load_sbml_cobra(sbml_file: Union[str, Path]) -> cobra.Model:
    """
    Loads a Pagerank-based Model from a COBRA SBML file.
    :param sbml_file: Path to json SBML file
    :type sbml_file: Union[str, Path]
    :return: Model object created from the COBRA model
    :rtype: Model
    """
    cobra_model = cobra.io.read_sbml_model(str(sbml_file))
    return cobra_model


def load_mat_cobra(mat_file: Union[str, Path]) -> cobra.Model:
    """
    Loads a Pagerank-based Model from a COBRA SBML file.
    :param mat_file: Path to json SBML file
    :type mat_file: Union[str, Path]
    :return: COBRA model object
    :rtype: Model
    """
    cobra_model = cobra.io.load_matlab_model(str(mat_file))
    return cobra_model


def load_csv(
    csv_file: Union[Path, str], sep=",", reversibilities: Optional[list[bool]] = None
):
    """
    Load models from CSV file (uses Pandas).
    :param csv_file: Path to CSV file
    :type csv_file: Union[Path, str]
    :param sep: Column separator, defaults to ','
    :type sep: str, optional
    :param reversibilities: List of reversibilities, defaults to None
    :type reversibilities: List[bool], optional
    :return: Model of the CSV file
    :rtype: Model
    """
    if not isinstance(csv_file, Path):
        csv_file = Path(csv_file)
    dataframe = pd.read_csv(csv_file, sep=sep, index_col=0)
    stoich_matrix = dataframe.values
    met = list(dataframe.index)
    rxn = list(dataframe.columns)
    if not reversibilities:
        reversibilities = [False] * len(rxn)
    return Model(stoich_matrix, met, reversibilities)


def pickle_model(
    model: Model, file_path: Union[Path, str], pickle_args: Optional[dict] = None
) -> Path:
    """
    Write a gemcat model to a pickle file
    :param model: Model to pickle
    :type model: Model
    :param file_path: Path/name to save pickle file to
    :type file_path: Union[Path, str]
    :param pickle_args: Arguments to pass on to pickle.dump, defaults to None
    :type pickle_args: Optional[dict], optional
    :return: Path to created pickle file
    :rtype: Path
    """
    if pickle_args is None:
        pickle_args = {}
    pickle.dump(model, str(file_path), **pickle_args)
    return Path(file_path)


def load_pickled(
    file_path: Union[Path, str], pickle_args: Optional[dict] = None
) -> Model:
    """
    Load a pickled gemcat model
    :param file_path: Path to pickle file
    :type file_path: str
    :param pickle_args: Arguments to pass to pickle.load, defaults to None
    :type pickle_args: Optional[dict], optional
    :return: Loaded model object
    :rtype: Model
    """
    if pickle_args is None:
        pickle_args = {}
    return pickle.load(str(file_path), **pickle_args)
