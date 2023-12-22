from pathlib import Path

import numpy as np
from fixtures import S_models, model_files_json, model_files_sbml, models

from gemcat import io


def test_convert_cobra_model(models, S_models):
    cobra_model = models["mini"]
    model = io.convert_cobra_model(cobra_model)
    assert np.allclose(model.stoichiometric_matrix, S_models["mini"])


def test_json_io(model_files_json, S_models):
    model_file = model_files_json["mini"]
    model = io.load_json_cobra(model_file)
    assert np.allclose(model.stoichiometric_matrix, S_models["mini"])


def test_sbml_io(model_files_sbml, S_models):
    model_file = model_files_sbml["mini"]
    model = io.load_sbml_cobra(model_file)
    assert np.allclose(model.stoichiometric_matrix, S_models["mini"])


def test_load_csv():
    filepath = Path("./tests/test_models/stoichiometry_trp.csv")
    model = io.load_csv(filepath, sep="\t")
    mets = [
        "M_Trp_c",
        "M_LKyn_c",
        "M_Serotonin",
        "M_3HAA",
        "M_3HKyn",
        "M_Trypta",
        "M_5HTrp",
        "M_FKyn",
        "M_Acms",
        "M_5HFKyn",
        "M_NaMN",
        "M_MTrypta",
        "M_FAA",
        "M_Quin",
    ]
    rxn_first = "(IL4I1_kcat_set_to_1)"
    rxn_last = "(AADAT_LKyn_kat4)"
    rxn_third = "(DDC_Trypta)"
    assert model.metabolite_names == mets
