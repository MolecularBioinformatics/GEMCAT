from pathlib import Path
from typing import Any

import requests

MODEL_URLS = {
    "recon3d": ("http://bigg.ucsd.edu/api/v2/models/Recon3D/download", "json"),
    "ratgem": (
        "https://github.com/SysBioChalmers/Rat-GEM/raw/refs/heads/main/model/Rat-GEM.mat",
        "mat",
    ),
}
SUPPORTED_MODELS = MODEL_URLS.keys()


def ratgem_processing(response: requests.Response) -> Any:
    """
    Parse request containing Rat-GEM .mat file.
    :param response: Response containing Rat-GEM model .mat binary
    :type response: requests.Response
    :return: Binary (.mat) representation of Rat-GEM
    :rtype: Any | Bytes
    """
    return response.content


def recon_processing(response: requests.Response) -> str:
    """
    Parse request containing Recon3D JSON
    :param response: Response containing Recon3D JSON
    :type response: requests.Response
    :return: String representation of Recon3D JSON
    :rtype: str
    """
    return response.text.replace("_AT", ".")


processing = {
    "recon3d": recon_processing,
    "ratgem": ratgem_processing,
}


# TODO: make singleton
class ModelManager:
    """
    ModelManager class in charge of management of auto-downloading, storing and retrieving model files
    """

    def __init__(self):
        """
        Initialize ModelManager
        """
        self.model_files = {
            name: f"{name}.{MODEL_URLS[name][1]}" for name in MODEL_URLS.keys()
        }
        self.model_path = Path("./models")
        self.model_file_paths = {
            name: self.model_path / model_file
            for name, model_file in self.model_files.items()
        }
        self.model_path.mkdir(exist_ok=True)
        self.managed_models_str = ", ".join(SUPPORTED_MODELS)

    def get_model(self, model: str) -> Path:
        """
        Retrieve a given model by its name
        :param model: Name given to the model (see allowed models above)
        :type model: str
        :return: Path to model file
        :rtype: Path
        """
        if not self.model_file_paths[model].exists():
            self.download_model(model)
        return self.model_file_paths[model]

    def download_model(self, model: str):
        """
        Download a given model by its name from a hardcoded URL
        :param model: Name given to the model (see allowed models above)
        :type model: str
        :raises ValueError: Raised if the model name is unknown
        :raises requests.HTTPError: Raised if the download fails
        """
        if not model in MODEL_URLS.keys():
            raise ValueError(
                f"Illegal model: {model} . Supported models are: {self.managed_models_str}"
            )
        try:
            response = requests.get(MODEL_URLS[model][0])
            response.raise_for_status()
        except requests.HTTPError as error:
            raise requests.HTTPError(f"Failed to download {model}") from error

        content = processing[model](response)
        save_file = self.model_file_paths[model]
        mode = "wb" if save_file.suffix == ".mat" else "w"
        with open(save_file, mode) as f:
            f.write(content)

    def get_managed_models(self) -> list[str]:
        """
        Return a list of the names of supported models
        :return: List of supported models' names
        :rtype: list[str]
        """
        return SUPPORTED_MODELS

    def get_managed_models_str(self) -> str:
        """
        Return all supported models in a string representation
        :return: String of supported models
        :rtype: str
        """
        return self.managed_models_str

    def wipe(self):
        """
        Wipes all previously downloaded model files and deletes the models folder
        :raises OSError: Raised if any file or the model folder cannot be deleted
        """
        try:
            for model in self.model_file_paths.values():
                if not model.exists():
                    continue
                model.unlink()
            self.model_path.rmdir()
        except Exception as error:
            raise OSError(
                "Could not delete model files, please delete manually in 'models' folder."
            ) from error
