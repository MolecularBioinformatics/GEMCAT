import cobra
from pathlib import Path
from typing import Union
from pyreporter import utils
from pyreporter.Model import Model
import pandas as pd
import pickle


def convert_cobra_model(cobra_model: cobra.Model) -> Model:
    """
    Converts a COBRA model into a PageRank-based Model.
    :param cobra_model: COBRA model object
    :type cobra_model: cobra.Model
    :return: Page-Rank based Model
    :rtype: Model
    """
    if not isinstance(cobra_model, cobra.Model):
        raise TypeError("COBRApy Model object expected")
    S = utils._get_stoich_matrix(cobra_model)
    rev = utils._get_reversibilities(cobra_model)
    metabolites = utils._get_metabolite_ids(cobra_model)
    model = Model(S, metabolites, rev)

    return model


def load_json_cobra(json_file: Union[str, Path]) -> Model:
    """
    Loads a PageRank-based Model from a COBRA json file.
    :param json_file: Path to json model file
    :type json_file: Union[str, Path]
    :return: Model object created from the COBRA model
    :rtype: Model
    """
    if isinstance(json_file, Path):
        json_file = json_file.as_posix()
    cobra_model = cobra.io.load_json_model(json_file)

    return convert_cobra_model(cobra_model)


def load_sbml_cobra(sbml_file: Union[str, Path]):
    """
    Loads a PageRank-based Model from a COBRA SBML file.
    :param sbml_file: Path to json SBML file
    :type sbml_file: Union[str, Path]
    :return: Model object created from the COBRA model
    :rtype: Model
    """
    if isinstance(sbml_file, Path):
        sbml_file = sbml_file.as_posix()
    cobra_model = cobra.io.read_sbml_model(sbml_file)

    return convert_cobra_model(cobra_model)


def load_csv(csv_file: Union[Path, str], sep=",", reversibilities=None):
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
    df = pd.read_csv(csv_file, sep=sep, index_col=0)
    S = df.values
    met = list(df.index)
    rxn = list(df.columns)
    if not reversibilities:
        reversibilities = [False] * len(rxn)
    S = Model(S, met, reversibilities)
    return S


def pickle_model(model: Model, file_name: str, pickle_args={}) -> Path:
    pickle.dump(model, file_name, **pickle_args)
    return Path(file_name)


def load_pickled(file_name: str, pickle_args={}) -> Model:
    return pickle.load(file_name, **pickle_args)
