from typing import Tuple

import numpy as np
import pyBigWig as bw

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def extract_bw(
    bw_path: str,
    chrom: str,
) -> Tuple[np.ndarray, int]:
    bwfile = bw.open(bw_path)
    if chrom not in bwfile.chroms():
        return None, None

    chrom_length = bwfile.chroms(chrom)
    data = bwfile.values(chrom, 0, chrom_length)
    return data, chrom_length
