import os
import shutil
from typing import Union

import numpy as np
import tiledb

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def create_tiledb_array_chrm(
    tiledb_uri_path: str,
    x_dim_length: int = None,
    y_dim_length: int = None,
    x_dim_name: str = "base",
    y_dim_name: str = "sample_index",
    matrix_attr_name: str = "data",
    x_dim_dtype: np.dtype = np.uint32,
    y_dim_dtype: np.dtype = np.uint32,
    matrix_dim_dtype: np.dtype = np.uint32,
    is_sparse: bool = False,
):
    """Create a TileDB file with the provided attributes to persistent storage.

    This will materialize the array directory and all
    related schema files.

    Args:
        tiledb_uri_path:
            Path to create the array TileDB file.

        x_dim_length:
            Number of entries along the x/fastest-changing dimension.
            e.g. Number of cells.
            Defaults to None, in which case, the max integer value of
            ``x_dim_dtype`` is used.

        y_dim_length:
            Number of entries along the y dimension.
            e.g. Number of genes.
            Defaults to None, in which case, the max integer value of
            ``y_dim_dtype`` is used.

        x_dim_name:
            Name for the x-dimension.
            Defaults to "cell_index".

        y_dim_name:
            Name for the y-dimension.
            Defaults to "gene_index".

        matrix_attr_name:
            Name for the attribute in the array.
            Defaults to "data".

        x_dim_dtype:
            NumPy dtype for the x-dimension.
            Defaults to np.uint32.

        y_dim_dtype:
            NumPy dtype for the y-dimension.
            Defaults to np.uint32.

        matrix_dim_dtype:
            NumPy dtype for the values in the matrix.
            Defaults to np.uint32.

        is_sparse:
            Whether the matrix is sparse.
            Defaults to False.
    """

    if x_dim_length is None:
        x_dim_length = np.iinfo(x_dim_dtype).max

    if y_dim_length is None:
        y_dim_length = np.iinfo(y_dim_dtype).max

    xdim = tiledb.Dim(name=x_dim_name, domain=(0, x_dim_length - 1), dtype=x_dim_dtype)
    ydim = tiledb.Dim(name=y_dim_name, domain=(0, y_dim_length - 1), dtype=y_dim_dtype)

    dom = tiledb.Domain(xdim, ydim)

    # expecting counts
    tdb_attr = tiledb.Attr(
        name=matrix_attr_name,
        dtype=matrix_dim_dtype,
        filters=tiledb.FilterList([tiledb.GzipFilter()]),
    )

    schema = tiledb.ArraySchema(domain=dom, sparse=is_sparse, attrs=[tdb_attr])

    if os.path.exists(tiledb_uri_path):
        shutil.rmtree(tiledb_uri_path)

    tiledb.Array.create(tiledb_uri_path, schema)

    tdbfile = tiledb.open(tiledb_uri_path, "w")
    tdbfile.close()


def write_array_to_tiledb(
    tiledb_array_uri: Union[str, tiledb.SparseArray],
    data: np.ndarray,
    x_idx: np.ndarray,
    y_idx: int,
    value_dtype: np.dtype = np.uint32,
):
    tiledb_fp = tiledb_array_uri
    if isinstance(tiledb_array_uri, str):
        tiledb_fp = tiledb.open(tiledb_array_uri, "w")

    if not isinstance(data, (np.ndarray)):
        raise TypeError("'data' is not an `ndarray`.")

    tiledb_fp[:x_idx, y_idx] = data.astype(value_dtype)


def optimize_tiledb_array(tiledb_array_uri: str, verbose: bool = True):
    """Consolidate TileDB fragments."""
    if verbose:
        print(f"Optimizing {tiledb_array_uri}")

    frags = tiledb.array_fragments(tiledb_array_uri)
    if verbose:
        print("Fragments before consolidation: {}".format(len(frags)))

    cfg = tiledb.Config()
    cfg["sm.consolidation.step_min_frags"] = 1
    cfg["sm.consolidation.step_max_frags"] = 200
    tiledb.consolidate(tiledb_array_uri, config=cfg)
    tiledb.vacuum(tiledb_array_uri)

    frags = tiledb.array_fragments(tiledb_array_uri)
    if verbose:
        print("Fragments after consolidation: {}".format(len(frags)))
