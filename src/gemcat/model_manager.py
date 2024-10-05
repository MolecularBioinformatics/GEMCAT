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
    return response.content


def recon_processing(response: requests.Response) -> str:
    return response.text.replace("_AT", ".")


processing = {
    "recon3d": recon_processing,
    "ratgem": ratgem_processing,
}


# TODO: make singleton
class ModelManager:
    def __init__(self):
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
        if not self.model_file_paths[model].exists():
            self.download_model(model)
        return self.model_file_paths[model]

    def download_model(self, model: str):
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

    def get_managed_models(self):
        return SUPPORTED_MODELS

    def get_managed_models_str(self):
        return self.managed_models_str

    def wipe(self):
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
