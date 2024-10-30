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


def extract_bw_intervals(
    bw_path: str,
    chrom: str,
) -> Optional[pd.DataFrame]:
    bwfile = bw.open(bw_path)
    if chrom not in bwfile.chroms():
        return None

    data = pd.DataFrame(bwfile.intervals(chrom), columns=["start", "end", "value"])
    data["chrom"] = chrom
    data_nnz = data[data["value"] > 0]

    data_nnz["tmp_group"] = (data_nnz["value"] != data_nnz["value"].shift()).cumsum()
    result = data_nnz.groupby("tmp_group").agg(
        {"start": "first", "end": "last", "value": "first"}
    )
    return result
