from typing import Optional, Tuple

import numpy as np
import pandas as pd
import pyBigWig as bw

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def extract_bw_values(
    bw_path: str,
    chrom: str,
) -> Tuple[np.ndarray, int]:
    bwfile = bw.open(bw_path)
    if chrom not in bwfile.chroms():
        return None, None

    chrom_length = bwfile.chroms(chrom)
    data = bwfile.values(chrom, 0, chrom_length)
    return np.array(data), chrom_length


# def extract_aspd():
#     df_data = pd.DataFrame(data, columns=["start", "end", "value"])
#     df_data["chrom"] = row.chrom
#     data_nnz = df_data[df_data["value"] > 0]

#     data_nnz["tmp_group"] = (
#         data_nnz["value"] != data_nnz["value"].shift()
#     ).cumsum()
#     results.append(
#         data_nnz.groupby("tmp_group").agg(
#             {"start": "first", "end": "last", "value": "first"}
#         )
#     )


def wrapper_extract_bw_values(
    bw_path: str,
    intervals: pd.DataFrame,
    agg_func: Optional[callable],
    val_dtype: np.dtype = np.float32,
    total_length: int = None,
) -> np.ndarray:
    if total_length is None:
        total_length = len(intervals)

    if agg_func is not None:
        return extract_bw_values_as_vec(bw_path=bw_path, intervals=intervals, agg_func=agg_func, val_dtype=val_dtype)
    else:
        return extract_bw_intervals_as_vec(
            bw_path=bw_path, intervals=intervals, val_dtype=val_dtype, total_length=total_length
        )


def extract_bw_values_as_vec(
    bw_path: str,
    intervals: pd.DataFrame,
    agg_func: Optional[callable] = None,
    val_dtype: np.dtype = np.float32,
) -> np.ndarray:
    """Extract data from BigWig for a given region and apply the aggregate function.

    Args:
        bw_path:
            Path to the BigWig file.

        intervals:
            List of intervals to extract.

        agg_func:
            Aggregate function to apply.
            Defaults to None.

        val_dtype:
            Dtype of the resulting array.

    Returns:
        A vector with length as the number of intervals,
        a value if the file contains the data for the corresponding
        region or ``np.nan`` if the region is not measured.
    """
    bwfile = bw.open(bw_path)

    results = []
    for row in intervals.itertuples():
        if row.seqnames in bwfile.chroms():
            try:
                data = bwfile.values(row.seqnames, row.starts, row.ends, numpy=True)
                if data is not None and len(data) != 0:
                    results.append(agg_func(data))
                else:
                    results.append(np.nan)
            except Exception as _:
                results.append(np.nan)
        else:
            results.append(np.nan)

    return np.array(results, dtype=val_dtype)


def _get_empty_array(size, val_dtype):
    out_array = np.empty(size, dtype=val_dtype)
    out_array.fill(np.nan)

    return out_array


def extract_bw_intervals_as_vec(
    bw_path: str,
    intervals: pd.DataFrame,
    total_length: int,
    val_dtype: np.dtype = np.float32,
) -> np.ndarray:
    """Extract data from BigWig for a given region.

    Args:
        bw_path:
            Path to the BigWig file.

        intervals:
            List of intervals to extract.

        total_length:
            Size of all the regions.

        val_dtype:
            Dtype of the resulting array.

    Returns:
        A vector with length as the number of intervals,
        a value if the file contains the data for the corresponding
        region or ``np.nan`` if the region is not measured.
    """
    bwfile = bw.open(bw_path)

    out_array = _get_empty_array(total_length, val_dtype=val_dtype)

    for row in intervals.itertuples():
        if row.seqnames in bwfile.chroms():
            try:
                data = bwfile.intervals(row.seqnames, row.starts, row.ends)
                tmp_out = _get_empty_array(row.ends - row.starts, val_dtype=val_dtype)
                if data is not None and len(data) != 0:
                    for d in data:
                        _strt = max(0, d[0] - row.starts)
                        _end = min(d[1] - row.starts, row.ends)
                        tmp_out[_strt:_end] = d[2]

                out_array[row.genarr_feature_index_start : row.genarr_feature_index_start + row.ends - row.starts] = (
                    tmp_out
                )
            except Exception as _:
                pass

    return out_array
