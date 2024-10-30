from dataclasses import dataclass
from typing import Dict

import numpy as np

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


@dataclass
class MatrixOptions:
    """Optional arguments for the ``matrix`` store for
    :py:func:`~genomicarrays.build_genomicarray.build_genomicarray`.

    Attributes:
        matrix_attr_name:
            Name of the matrix to be stored in the TileDB file.
            Defaults to "data".

        skip:
            Whether to skip generating matrix TileDB.
            Defaults to False.

        dtype:
            NumPy dtype for the values in the matrix.
            Defaults to np.uint16.

            Note: make sure the matrix values fit
            within the range limits of chosen-dtype.
    """

    skip: bool = False
    matrix_attr_name: str = "data"
    dtype: np.dtype = np.uint16


@dataclass
class SampleMetadataOptions:
    """Optional arguments for the ``sample`` store for
    :py:func:`~genomicarrays.build_genomicarray.build_genomicarray`.

    Attributes:
        skip:
            Whether to skip generating sample TileDB.
            Defaults to False.

        dtype:
            NumPy dtype for the sample dimension.
            Defaults to np.uint32.

            Note: make sure the number of samples fit
            within the integer limits of chosen dtype.

        tiledb_store_name:
            Name of the TileDB file.
            Defaults to "sample_metadata".

        column_types:
            A dictionary containing column names as keys
            and the value representing the type to in
            the TileDB.

            If `None`, all columns are cast as 'ascii'.
    """

    skip: bool = False
    dtype: np.dtype = np.uint32
    tiledb_store_name: str = "sample_metadata"
    column_types: Dict[str, np.dtype] = None
