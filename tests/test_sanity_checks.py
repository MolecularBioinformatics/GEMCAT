import numpy as np
import pandas as pd
from fixtures import models

import gemcat as pr


def test_sanity_double_all(models):
    genes = [f"G{i}" for i in range(1, 6)]
    genes_mapped_single = pd.Series(dict(zip(genes, [1.0] * 5)))
    genes_mapped_double = pd.Series(dict(zip(genes, [2.0] * 5)))
    mini = models["mini"]
    result_single = pr.workflows.workflow_avg_single(mini, genes_mapped_single)
    result_double = pr.workflows.workflow_avg_single(mini, genes_mapped_double)
    assert np.allclose(result_single.values, result_double.values)
