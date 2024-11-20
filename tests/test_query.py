import tempfile

import anndata
import numpy as np
import pandas as pd
import pytest
import tiledb
from genomicarrays import (
    GenomicArrayDataset,
    build_genomicarray,
    FeatureAnnotationOptions,
)

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def test_query_agg():
    tempdir = tempfile.mkdtemp()

    strts = np.arange(300, 600, 20)
    features = pd.DataFrame(
        {"seqnames": ["chr1"] * 15, "starts": strts, "ends": strts + 15}
    )

    dataset = build_genomicarray(
        output_path=tempdir,
        files=["tests/data/test1.bw", "tests/data/test2.bw"],
        features=features,
        genome_fasta="tests/data/test.fa",
        feature_annotation_options=FeatureAnnotationOptions(
            aggregate_function=np.nanmean
        ),
    )

    assert dataset is not None
    assert isinstance(dataset, GenomicArrayDataset)

    cd = GenomicArrayDataset(dataset_path=tempdir)

    ffp = tiledb.open(f"{tempdir}/feature_annotation", "r")

    features_rt = ffp.df[:]
    assert np.allclose(np.array(features["starts"]), np.array(features_rt["starts"]))
    assert np.allclose(np.array(features["ends"]), np.array(features_rt["ends"]))
    assert all([a == b for a, b in zip(features["seqnames"], features_rt["seqnames"])])
    assert "genarr_feature_index" in features_rt.columns

    assert cd.shape == (15, 2)
    assert len(cd) == 15

    result1 = cd[:, 0]

    assert result1 is not None
    assert result1.matrix.shape == (len(features), 1)

    assert np.allclose(
        result1.matrix.flatten(),
        np.repeat(1, 15),
    )

    result2 = cd[:, 1]
    assert result2.matrix.shape == (len(features), 1)

    expected = np.concatenate((np.repeat(0.5, 10), np.repeat(np.nan, 5)))

    assert np.allclose(
        result2.matrix.flatten()[:10],
        expected[:10],
    )

    assert cd.get_sample_metadata_columns() == [
        "genarr_sample",
    ]
    assert len(cd.get_sample_metadata_column("genarr_sample")) == 2
    assert len(cd.get_sample_subset("genarr_sample == 'sample_1'")) == 1

    assert sorted(cd.get_feature_annotation_columns()) == sorted(
        [
            "ends",
            "genarr_feature_index",
            "genarr_feature_start_index",
            "genarr_feature_end_index",
            "seqnames",
            "sequences",
            "starts",
            "widths",
        ]
    )
    assert len(cd.get_feature_annotation_column("genarr_feature_index")) == 15
    assert len(cd.get_feature_subset("genarr_feature_index == 1")) == 1

    result1 = cd.get_slice(slice(0, 5), slice(None))
    assert result1 is not None
    assert result1.matrix.shape == (6, 2)

    assert result1.to_anndata() is not None
    assert result1.to_rangedsummarizedexperiment() is not None


def test_query_noagg():
    tempdir = tempfile.mkdtemp()

    strts = np.arange(300, 600, 20)
    features = pd.DataFrame(
        {"seqnames": ["chr1"] * 15, "starts": strts, "ends": strts + 15}
    )

    dataset = build_genomicarray(
        output_path=tempdir,
        files=["tests/data/test1.bw", "tests/data/test2.bw"],
        features=features,
        genome_fasta="tests/data/test.fa",
    )

    assert dataset is not None
    assert isinstance(dataset, GenomicArrayDataset)

    cd = GenomicArrayDataset(dataset_path=tempdir)

    ffp = tiledb.open(f"{tempdir}/feature_annotation", "r")

    features_rt = ffp.df[:]
    assert np.allclose(np.array(features["starts"]), np.array(features_rt["starts"]))
    assert np.allclose(np.array(features["ends"]), np.array(features_rt["ends"]))
    assert all([a == b for a, b in zip(features["seqnames"], features_rt["seqnames"])])
    assert "genarr_feature_index" in features_rt.columns

    assert cd.shape == (15, 2)
    assert len(cd) == 15

    result1 = cd[:, 0]

    assert result1 is not None
    assert result1.matrix.shape == (225, 1)

    assert np.all(np.isnan(result1.matrix.flatten()))

    result2 = cd[:, 1]
    assert result2.matrix.shape == (225, 1)

    assert np.all(np.isnan(result2.matrix.flatten()))

    assert cd.get_sample_metadata_columns() == [
        "genarr_sample",
    ]
    assert len(cd.get_sample_metadata_column("genarr_sample")) == 2
    assert len(cd.get_sample_subset("genarr_sample == 'sample_1'")) == 1

    assert sorted(cd.get_feature_annotation_columns()) == sorted(
        [
            "ends",
            "genarr_feature_index",
            "genarr_feature_start_index",
            "genarr_feature_end_index",
            "seqnames",
            "sequences",
            "starts",
            "widths",
        ]
    )
    assert len(cd.get_feature_annotation_column("genarr_feature_index")) == 15
    assert len(cd.get_feature_subset("genarr_feature_index == 1")) == 1

    result1 = cd.get_slice(slice(0, 5), slice(None))
    assert result1 is not None
    assert result1.matrix.shape == (90, 2)

    assert result1.to_anndata() is not None
    assert result1.to_rangedsummarizedexperiment() is not None
