"""Build the `GenomicArrayDatset`.

The `GenomicArrayDatset` method is designed to store genomic range-based 
datasets from BigWig, BigBed or other similar files.

Example:

    .. code-block:: python

        import pyBigWig as bw
        import numpy as np
        import tempfile
        from genomicarrays import build_genomicarray, MatrixOptions

        # Create a temporary directory
        tempdir = tempfile.mkdtemp()

        # Read BigWig objects
        bw1 = bw.open("path/to/object1.bw", "r")
        # or just provide the path
        bw2 = "path/to/object2.bw"

        # Build GenomicArray
        dataset = build_genomicarray(
            output_path=tempdir,
            files=[bw1, bw2],
            matrix_options=MatrixOptions(dtype=np.float32),
        )
"""

import os
import warnings
from typing import Dict, Union

import pyBigWig as bw
import numpy as np
import pandas as pd
import tiledb
from multiprocessing import Pool

from . import build_options as bopt
from . import buildutils_tiledb_array as uta
from cellarr import buildutils_tiledb_frame as utf
from . import utils_bw as ubw
# from .GenomicArrayDataset import GenomicArrayDataset

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


# TODO: Accept files as a dictionary with names to each dataset.
# TODO: options to only extract specific regions.
def build_genomicarray(
    files: list,
    output_path: str,
    genome: Union[str, pd.DataFrame] = "hg38",
    # filter_to_regions: Union[str, pd.DataFrame],
    sample_metadata: Union[pd.DataFrame, str] = None,
    sample_metadata_options: bopt.SampleMetadataOptions = bopt.SampleMetadataOptions(),
    matrix_options: bopt.MatrixOptions = bopt.MatrixOptions(),
    optimize_tiledb: bool = True,
    num_threads: int = 1,
):
    """Generate the `GenomicArrayDatset`.

    All files are expected to be consistent and any modifications
    to make them consistent is outside the scope of this function
    and package.

    Args:
        files:
            List of file paths to `BigWig` files.

        output_path:
            Path to where the output TileDB files should be stored.

        sample_metadata:
            A :py:class:`~pandas.DataFrame` containing the sample
            metadata for each file in ``files``. Hences the number of rows
            in the dataframe must match the number of ``files``.

            Alternatively, may provide path to the file containing a
            concatenated sample metadata across all BigWig files. In this case,
            the first row is expected to contain the column names.

            Additionally, the order of rows is expected to be in the same
            order as the input list of ``files``.

            Defaults to `None`, in which case, we create a simple sample
            metadata dataframe containing the list of datasets, aka 
            each BigWig files. Each dataset is named as ``sample_{i}``
            where `i` refers to the index position of the object in ``files``.

        sample_metadata_options:
            Optional parameters when generating ``sample_metadata`` store.

        matrix_options:
            Optional parameters when generating ``matrix`` store.

        optimize_tiledb:
            Whether to run TileDB's vaccum and consolidation (may take long).

        num_threads:
            Number of threads.
            Defaults to 1.
    """
    if not os.path.isdir(output_path):
        raise ValueError("'output_path' must be a directory.")
    
    ####
    ## Writing the sample metadata file
    ####

    if isinstance(genome, str):
        chrom_sizes = pd.read_csv("https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes", sep="\t", header=None, names=["chrom", "length"])
    elif isinstance(genome, pd.DataFrame):
        chrom_sizes = genome
    else:
        raise TypeError("'genome' is not an expected type (either 'str' or 'Dataframe').")

    ####
    ## Writing the sample metadata file
    ####
    _samples = []
    for idx, _ in enumerate(files):
        _samples.append(f"sample_{idx + 1}")

    if sample_metadata is None:
        warnings.warn(
            "Sample metadata is not provided, each dataset in 'files' is considered a sample",
            UserWarning,
        )

        sample_metadata = pd.DataFrame({"genarr_sample": _samples})
    elif isinstance(sample_metadata, str):
        sample_metadata = pd.read_csv(sample_metadata, header=0)
        if "genarr_sample" not in sample_metadata.columns:
            sample_metadata["genarr_sample"] = _samples
    elif isinstance(sample_metadata, pd.DataFrame):
        if "genarr_sample" not in sample_metadata.columns:
            sample_metadata["genarr_sample"] = _samples
    else:
        raise TypeError("'sample_metadata' is not an expected type.")

    if not sample_metadata_options.skip:
        _col_types = utf.infer_column_types(
            sample_metadata, sample_metadata_options.column_types
        )

        _sample_output_uri = (
            f"{output_path}/{sample_metadata_options.tiledb_store_name}"
        )
        utf.create_tiledb_frame_from_dataframe(
            _sample_output_uri, sample_metadata, column_types=_col_types
        )

        if optimize_tiledb:
            uta.optimize_tiledb_array(_sample_output_uri)

    ####
    ## Writing the genomic ranges file
    ####
    if not matrix_options.skip:
        tiledb.group_create(f"{output_path}/ranges")

        _chrm_group_base = f"{output_path}/ranges"
        for idx, seq in chrom_sizes.iterrows():
            uta.create_tiledb_array(
                f"{_chrm_group_base}/{seq['chrom']}",
                x_dim_length = seq['length'],
                y_dim_length = len(files),
                matrix_attr_name = matrix_options.matrix_attr_name,
                matrix_dim_dtype = matrix_options.dtype,
            )

        if optimize_tiledb:
            uta.optimize_tiledb_array(_counts_uri)

    # return GenomicArrayDataset(
    #     dataset_path=output_path,
    #     sample_metadata_uri=sample_metadata_options.tiledb_store_name,
    #     cell_metadata_uri=cell_metadata_options.tiledb_store_name,
    #     gene_annotation_uri=gene_annotation_options.tiledb_store_name,
    #     matrix_tdb_uri=matrix_options.tiledb_store_name,
    # )

def _range_writer(outpath, bwpath, bwidx, chrm):


def _wrapper_extract_info(args):
    outpath, bwpath, bwidx, chrm = args
    return _range_writer(outpath, bwpath, bwidx, chrm)


def extract_anndata_info(
    h5ad_or_adata: List[Union[str, anndata.AnnData]],
    var_feature_column: str = "index",
    obs_subset_columns: dict = None,
    num_threads: int = 1,
):
    with Pool(num_threads) as p:
        _args = [
            (file_info, var_feature_column, obs_subset_columns)
            for file_info in h5ad_or_adata
        ]
        return p.map(_wrapper_extract_info, _args)