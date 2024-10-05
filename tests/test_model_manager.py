from pathlib import Path

from pytest import fixture

from gemcat.model_manager import ModelManager

model_path = Path("./models")


def test_model_manager_setup():
    try:
        mm = ModelManager()
        assert isinstance(mm, ModelManager)
    except Exception:
        raise AssertionError()


def test_wipe():
    try:
        mm = ModelManager()
        mm.model_path.mkdir(exist_ok=True)
        fake_model = mm.model_path / "recon3d.json"
        with open(fake_model, "w") as f:
            f.write("delete me")
        assert (mm.model_path / "recon3d.json").is_file()
        mm.wipe()
        assert not (mm.model_path / "recon3d.json").is_file()
    finally:
        if fake_model.exists():
            fake_model.unlink()


@fixture
def ensure_empty_models():
    mm = ModelManager()
    mm.wipe()


def test_model_manager_recon(ensure_empty_models):
    mm = ModelManager()
    mm.download_model("recon3d")
    assert (mm.model_path / "recon3d.json").exists()


def test_model_manager_ratgem(ensure_empty_models):
    mm = ModelManager()
    mm.download_model("ratgem")
    assert (mm.model_path / "ratgem.mat").exists()


def test_model_manager_get_download(ensure_empty_models):
    mm = ModelManager()
    output = mm.get_model("ratgem")
    assert output == mm.model_path / "ratgem.mat"
    assert (output).exists()


def test_model_manager_get_existing(ensure_empty_models):
    try:
        mm = ModelManager()
        output = mm.model_path / "recon3d.json"
        with open(output, "w") as f:
            f.write("I am a model.")
        result = mm.get_model("recon3d")
        assert result == output
        assert result.exists()
    finally:
        mm.wipe()
