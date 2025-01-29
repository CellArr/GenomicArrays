import tempfile

import numpy as np
import pandas as pd
import tiledb

from genomicarrays import build_genomicarray, FeatureAnnotationOptions
from agg_test import aggfunc

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def test_ingest_bigwigs():
    tempdir = tempfile.mkdtemp()

    strts = np.arange(300, 600, 20)
    features = pd.DataFrame(
        {"seqnames": ["chr1"] * 15, "starts": strts, "ends": strts + 15}
    )

    build_genomicarray(
        output_path=tempdir,
        files=["tests/data/test1.bw", "tests/data/test2.bw"],
        features=features,
        genome_fasta="tests/data/test.fa"
    )

    cfp = tiledb.open(f"{tempdir}/coverage", "r")
    ffp = tiledb.open(f"{tempdir}/feature_annotation", "r")

    features = ffp.df[:]

    assert len(features) == 15

    assert np.all(np.isnan(cfp.multi_index[:5, 0]["data"].flatten()))
    assert np.all(np.isnan(cfp.multi_index[:5, 1]["data"].flatten()))

    sfp = tiledb.open(f"{tempdir}/sample_metadata", "r")
    samples = sfp.df[:]
    assert len(samples) == 2

def test_ingest_bigwigs_agg():
    tempdir = tempfile.mkdtemp()

    strts = np.arange(300, 600, 20)
    features = pd.DataFrame(
        {"seqnames": ["chr1"] * 15, "starts": strts, "ends": strts + 15}
    )

    build_genomicarray(
        output_path=tempdir,
        files=["tests/data/test1.bw", "tests/data/test2.bw"],
        features=features,
        genome_fasta="tests/data/test.fa",
        feature_annotation_options=FeatureAnnotationOptions(aggregate_function=np.nanmean)
    )

    cfp = tiledb.open(f"{tempdir}/coverage", "r")
    ffp = tiledb.open(f"{tempdir}/feature_annotation", "r")

    features = ffp.df[:]

    assert len(features) == 15
    assert len(cfp) == 15

    assert np.allclose(
        cfp.multi_index[0:5, 0]["data"],
        np.repeat(1, 5),
    )
    assert np.allclose(
        cfp.multi_index[0:5, 1]["data"],
        np.repeat(0.5, 5),
    )

    sfp = tiledb.open(f"{tempdir}/sample_metadata", "r")
    samples = sfp.df[:]
    assert len(samples) == 2

def test_ingest_bigwigs_agg_not_scalar():
    tempdir = tempfile.mkdtemp()

    strts = np.arange(300, 600, 20)
    features = pd.DataFrame(
        {"seqnames": ["chr1"] * 15, "starts": strts, "ends": strts + 15}
    )

    build_genomicarray(
        output_path=tempdir,
        files=["tests/data/test1.bw", "tests/data/test2.bw"],
        features=features,
        genome_fasta="tests/data/test.fa",
        feature_annotation_options=FeatureAnnotationOptions(aggregate_function=aggfunc, expected_agg_function_length=2),
    )

    cfp = tiledb.open(f"{tempdir}/coverage", "r")
    ffp = tiledb.open(f"{tempdir}/feature_annotation", "r")

    features = ffp.df[:]

    assert len(features) == 15
    assert len(cfp) == 30

    assert np.allclose(
        cfp.multi_index[0:5, 0]["data"].flatten(),
        np.array([1,10,1,10,1,10]),
    )
    assert np.allclose(
        cfp.multi_index[0:5, 1]["data"].flatten(),
        np.array([0.5, 5,0.5, 5,0.5, 5]),
    )

    sfp = tiledb.open(f"{tempdir}/sample_metadata", "r")
    samples = sfp.df[:]
    assert len(samples) == 2