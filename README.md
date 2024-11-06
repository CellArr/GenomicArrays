<!-- These are examples of badges you might want to add to your README:
     please update the URLs accordingly

[![Built Status](https://api.cirrus-ci.com/github/<USER>/GenomicArrays.svg?branch=main)](https://cirrus-ci.com/github/<USER>/GenomicArrays)
[![ReadTheDocs](https://readthedocs.org/projects/GenomicArrays/badge/?version=latest)](https://GenomicArrays.readthedocs.io/en/stable/)
[![Coveralls](https://img.shields.io/coveralls/github/<USER>/GenomicArrays/main.svg)](https://coveralls.io/r/<USER>/GenomicArrays)
[![PyPI-Server](https://img.shields.io/pypi/v/GenomicArrays.svg)](https://pypi.org/project/GenomicArrays/)
[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/GenomicArrays.svg)](https://anaconda.org/conda-forge/GenomicArrays)
[![Monthly Downloads](https://pepy.tech/badge/GenomicArrays/month)](https://pepy.tech/project/GenomicArrays)
[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter)](https://twitter.com/GenomicArrays)
-->

[![Project generated with PyScaffold](https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold)](https://pyscaffold.org/)

# Genomic Arrays based on TileDB

GenomicArrays is a Python package for converting genomic data from BigWig format to TileDB arrays.

## Installation

Install the package from [PyPI](https://pypi.org/project/genomicarrays/)

```sh
pip install genomicarrays
```

## Quick Start

### Build a `GenomicArray`

Building a `GenomicArray` generates 3 TileDB files in the specified output directory:

- `feature_annotation`: A TileDB file containing input feature intervals.
- `sample_metadata`: A TileDB file containing sample metadata, each BigWig file is considered a sample.
- A matrix TileDB file named by the `layer_matrix_name` parameter. This allows the package
to store multiple different matrices, e.g. 'coverage', 'some_computed_statistic', for the same interval,
and sample metadata attributes.

The organization is inspired by the [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) data structure. The TileDB matrix file is stored in a **features X samples** orientation.

![`GenomicArray` structure](./assets/genarr.png "GenomicArray")

To build a `GenomicArray` from a collection of `BigWig` files:

```python
import numpy as np
import tempfile
import genomicarrays as garr

# Create a temporary directory, this is where the
# output files are created. Pick your location here.
tempdir = tempfile.mkdtemp()

# List BigWig paths
bw_dir = "your/biwig/dir"
files = os.listdir(bw_dir)
bw_files = [f"{bw_dir}/{f}" for f in files]

features = pd.DataFrame({
     "seqnames": ["chr1", "chr1"],
     "starts": [1000, 2000],
     "ends": [1500, 2500]
})

# Build GenomicArray
dataset = garr.build_genomicarray(
     files=bw_files,
     output_path=tempdir,
     features=features,
     # Specify a fasta file to extract sequences
     # for each region in features
     genome_fasta="path/to/genome.fasta",
     # agg function to summarize mutiple values
     # from bigwig within an input feature interval.
     feature_annotation_options=garr.FeatureAnnotationOptions(
        aggregate_function = np.nanmean
     ),
     # for parallel processing multiple bigwig files
     num_threads=4
)
```

The build process stores missing intervals from a bigwig file as `np.nan`. The
default is to choose an aggregate functions that works with `np.nan`.

### Query a `GenomicArrayDataset`

Users have the option to reuse the `dataset` object retuned when building the arrays or by creating a `GenomicArrayDataset` object by initializing it to the path where the files were created.

```python
# Create a GenomicArrayDataset object from the existing dataset
dataset = GenomicArrayDataset(dataset_path=tempdir)

# Query data for the first 10 regions across all samples
coverage_data = dataset[0:10, :]

print(expression_data.matrix)
print(expression_data.feature_annotation)
```

     ## output 1
     array([[1. , 0.5],
          [1. , 0.5],
          [1. , 0.5],
          [1. , 0.5],
          [1. , 0.5],
          [1. , 0.5],
          [1. , 0.5],
          [1. , 0.5],
          [1. , 0.5],
          [1. , 0.5],
          [1. , nan]], dtype=float32)

     ## output 2
     seqnames  starts  ends  genarr_feature_index
     0      chr1     300   315                     0
     1      chr1     320   335                     1
     2      chr1     340   355                     2
     3      chr1     360   375                     3
     4      chr1     380   395                     4
     5      chr1     400   415                     5
     6      chr1     420   435                     6
     7      chr1     440   455                     7
     8      chr1     460   475                     8
     9      chr1     480   495                     9
     10     chr1     500   515                    10


<!-- pyscaffold-notes -->

## Note

This project has been set up using PyScaffold 4.6. For details and usage
information on PyScaffold see https://pyscaffold.org/.
