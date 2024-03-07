from pathlib import Path

import pytest
import requests

model_path = Path("./tests/test_models/")
expression_path = Path("./tests/test_seq/")


def download_link(link: str, file_name: str, file_path: Path):
    """
    Download the file from a link to a given path and name.
    Skip with a print if it already exists at that location.
    :param link: Link to file.
    :type link: string
    :param file_name: file name
    :type file_name: string
    :param file_path: path to save file at
    :type file_path: Path
    """
    file_path = file_path / file_name
    if file_path.is_file():
        print(f"Skipping existing {file_name}")
        return
    response = requests.get(link)
    try:
        response.raise_for_status()
        with open(file_path, "w") as f:
            f.write(response.text)
    except requests.HTTPError as err:
        print(
            f"Failed to download {file_name} from {link}. Tests may fail as a result."
        )
        print(err)


def download_models():
    """
    Download all models used for testing
    """
    print("Downloading models for testing...")
    models_to_download = {
        "recon3d.xml": "http://bigg.ucsd.edu/static/models/Recon3D.xml",
        "recon3d.json": "http://bigg.ucsd.edu/static/models/Recon3D.json",
    }
    for file_name, link in models_to_download.items():
        download_link(link, file_name, model_path)
    print("Model downloads have finished.")


def download_expression():
    """
    Download all expression files used for testing
    """
    print("Downloading expression data for testing...")
    expression_to_download = {
        "test_prot_uc_vs_healthy.csv": "https://media.githubusercontent.com/media/MolecularBioinformatics/prm_manuscript/main/data/prot_uc_vs_healthy.csv",
        "test_rnaseq_HEK293_complete.csv": "https://media.githubusercontent.com/media/MolecularBioinformatics/prm_manuscript/main/data/rnaseq_HEK293_complete.csv",
        "test_rnaseq_slc25a51ko_vs_parental.csv": "https://media.githubusercontent.com/media/MolecularBioinformatics/prm_manuscript/main/data/rnaseq_slc25a51ko_vs_parental.csv",
    }
    for file_name, link in expression_to_download.items():
        download_link(link, file_name, expression_path)
    print("Expression data downloads have finished.")


def pytest_sessionstart():
    """
    Tasks to perform when a test session is started
    """
    download_models()
    download_expression()
