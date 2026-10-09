# Writes the regression fixture for GregorLueg/bixverse#252: one obs column per
# pandas column kind, as anndata >= 0.13 encodes them.
# uv run --python 3.12 --with anndata==0.13.4 --with scipy \
#   python data-raw/anndata_nullable_h5ad.py
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

n = 40
rng = np.random.default_rng(0)
obs = pd.DataFrame(
    {
        "donor": [f"donor{i}" for i in range(n)],
        "sample": pd.array([f"s{i}" for i in range(n)], dtype="string"),
        "batch": ["A", "B"] * (n // 2),
        "n_doublets": pd.array([1, None] * (n // 2), dtype="Int64"),
        "is_doublet": pd.array([True, None] * (n // 2), dtype="boolean"),
        "age": np.arange(n),
        "passed_qc": np.arange(n) % 2 == 0,
        "score": rng.random(n),
    },
    index=[f"cell{i}" for i in range(n)],
)
var = pd.DataFrame(
    {"symbol": [f"SYM{i}" for i in range(5)]},
    index=[f"gene{i}" for i in range(5)],
)
X = sp.csr_matrix(rng.poisson(5, (n, 5)).astype(np.float32))
ad.AnnData(X, obs=obs, var=var).write_h5ad(
    "inst/tinytest/synthetic_data/anndata_0_13_nullable.h5ad"
)
