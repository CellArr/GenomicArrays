"""Query the `GenomicArrayDataset`.

This class provides methods to access the directory containing the
generated TileDB files usually using the
:py:func:`~genomicarray.build_genomicarray.build_genomicarray`.

Example:

    .. code-block:: python

        from genomicarray import GenomicArrayDataset

        garr = GenomicArrayDataset(dataset_path="/path/to/genomicarray/dir")
        result1 = garr[0:10, 0]

        print(result1)
"""

import os
from typing import List, Sequence, Union

import pandas as pd
import tiledb

from . import queryutils_tiledb_frame as qtd
from .GenomicArrayDatasetSlice import GenomicArrayDatasetSlice

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


class GenomicArrayDataset:
    """A class that represent a collection of features and their associated coverage in a TileDB backed store."""

    def __init__(
        self,
        dataset_path: str,
        matrix_tdb_uri: str = "coverage",
        feature_annotation_uri: str = "feature_annotation",
        sample_metadata_uri: str = "sample_metadata",
    ):
        """Initialize a ``GenomicArrayDataset``.

        Args:
            dataset_path:
                Path to the directory containing the TileDB stores.
                Usually the ``output_path`` from the
                :py:func:`~genomicarray.build_genomicarray.build_genomicarray`.

            matrix_tdb_uri:
                Relative path to matrix store.

            feature_annotation_uri:
                Relative path to feature annotation store.

            sample_metadata_uri:
                Relative path to sample metadata store.
        """

        if not os.path.isdir(dataset_path):
            raise ValueError("'dataset_path' is not a directory.")

        self._dataset_path = dataset_path
        # TODO: Maybe switch to on-demand loading of these objects
        self._matrix_tdb_tdb = tiledb.open(f"{dataset_path}/{matrix_tdb_uri}", "r")
        self._feature_annotation_tdb = tiledb.open(f"{dataset_path}/{feature_annotation_uri}", "r")
        self._sample_metadata_tdb = tiledb.open(f"{dataset_path}/{sample_metadata_uri}", "r")

    def __del__(self):
        self._matrix_tdb_tdb.close()
        self._feature_annotation_tdb.close()
        self._sample_metadata_tdb.close()

    ####
    ## Subset methods for the `feature_annotation` TileDB file.
    ####
    def get_feature_annotation_columns(self) -> List[str]:
        """Get annotation column names from ``feature_annotation`` store.

        Returns:
            List of available annotations.
        """
        return qtd.get_schema_names_frame(self._feature_annotation_tdb)

    def get_feature_annotation_column(self, column_name: str) -> pd.DataFrame:
        """Access a column from the ``feature_annotation`` store.

        Args:
            column_name:
                Name of the column or attribute. Usually one of the column names
                from of :py:meth:`~get_feature_annotation_columns`.

        Returns:
            A list of values for this column.
        """
        res = qtd.get_a_column(self._feature_annotation_tdb, column_name=column_name)
        return res[column_name]

    def get_feature_annotation_index(self) -> List[str]:
        """Get index of the ``feature_annotation`` store.

        Returns:
            List of feature ids.
        """
        res = qtd.get_a_column(self._feature_annotation_tdb, "genarr_feature_index")
        return res["genarr_feature_index"].tolist()

    def get_feature_subset(self, subset: Union[slice, List[str], tiledb.QueryCondition], columns=None) -> pd.DataFrame:
        """Slice the ``feature_annotation`` store.

        Args:
            subset:
                A list of integer indices to subset the ``feature_annotation``
                store.

                Alternatively, may provide a
                :py:class:`tiledb.QueryCondition` to query the store.

                Alternatively, may provide a list of strings to match with
                the index of ``feature_annotation`` store.

            columns:
                List of specific column names to access.

                Defaults to None, in which case all columns are extracted.

        Returns:
            A pandas Dataframe of the subset.
        """

        if isinstance(columns, str):
            columns = [columns]

        if columns is None:
            columns = self.get_feature_annotation_columns()
        else:
            _not_avail = []
            for col in columns:
                if col not in self.get_feature_annotation_columns():
                    _not_avail.append(col)

            if len(_not_avail) > 0:
                raise ValueError(f"Columns '{', '.join(_not_avail)}' are not available.")

        if qtd._is_list_strings(subset):
            subset = self._get_indices_for_gene_list(subset)

        return qtd.subset_frame(self._feature_annotation_tdb, subset=subset, columns=columns)

    ####
    ## Subset methods for the `sample_metadata` TileDB file.
    ####
    def get_sample_metadata_columns(self) -> List[str]:
        """Get column names from ``sample_metadata`` store.

        Returns:
            List of available metadata columns.
        """
        return qtd.get_schema_names_frame(self._sample_metadata_tdb)

    def get_sample_metadata_column(self, column_name: str) -> pd.DataFrame:
        """Access a column from the ``sample_metadata`` store.

        Args:
            column_name:
                Name of the column or attribute. Usually one of the column names
                from of :py:meth:`~get_sample_metadata_columns`.

        Returns:
            A list of values for this column.
        """
        res = qtd.get_a_column(self._sample_metadata_tdb, column_name=column_name)
        return res[column_name]

    def get_sample_subset(self, subset: Union[slice, tiledb.QueryCondition], columns=None) -> pd.DataFrame:
        """Slice the ``sample_metadata`` store.

        Args:
            subset:
                A list of integer indices to subset the ``sample_metadata``
                store.

                Alternatively, may also provide a
                :py:class:`tiledb.QueryCondition` to query the store.

            columns:
                List of specific column names to access.

                Defaults to None, in which case all columns are extracted.

        Returns:
            A pandas Dataframe of the subset.
        """
        if isinstance(columns, str):
            columns = [columns]

        if columns is None:
            columns = self.get_sample_metadata_columns()
        else:
            _not_avail = []
            for col in columns:
                if col not in self.get_sample_metadata_columns():
                    _not_avail.append(col)

            if len(_not_avail) > 0:
                raise ValueError(f"Columns '{', '.join(_not_avail)}' are not available.")

        return qtd.subset_frame(self._sample_metadata_tdb, subset=subset, columns=columns)

    ####
    ## Subset methods for the `matrix` TileDB file.
    ####
    def get_matrix_subset(self, subset: Union[int, Sequence, tuple]) -> pd.DataFrame:
        """Slice the ``matrix`` store.

        Args:
            subset:
                Any `slice`supported by TileDB's array slicing.
                For more info refer to
                <TileDB docs https://docs.tiledb.com/main/how-to/arrays/reading-arrays/basic-reading>_.

        Returns:
            A pandas Dataframe of the subset.
        """
        if isinstance(subset, (str, int)):
            return qtd.subset_array(
                self._matrix_tdb_tdb,
                subset,
                slice(None),
                shape=(len(subset), self.shape[1]),
            )

        if isinstance(subset, tuple):
            if len(subset) == 0:
                raise ValueError("At least one slicing argument must be provided.")

            if len(subset) == 1:
                return qtd.subset_array(
                    self._matrix_tdb_tdb,
                    subset[0],
                    slice(None),
                    shape=(len(subset[0]), self.shape[1]),
                )
            elif len(subset) == 2:
                return qtd.subset_array(
                    self._matrix_tdb_tdb,
                    subset[0],
                    subset[1],
                    shape=(len(subset[0]), len(subset[1])),
                )
            else:
                raise ValueError(f"`{type(self).__name__}` only supports 2-dimensional slicing.")

    ####
    ## Subset methods by cell and gene dimensions.
    ####
    def get_slice(
        self,
        feature_subset: Union[slice, int],
        sample_subset: Union[slice, List[str], tiledb.QueryCondition],
    ) -> GenomicArrayDatasetSlice:
        """Subset a ``GenomicArrayDataset``.

        Args:
            sample_subset:
                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the columns (or samples) to retain.

            feature_subset:
                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the rows (or features/genes) to retain.

        Returns:
            A :py:class:`~genomicarray.GenomicArrayDatasetSlice.GenomicArrayDatasetSlice` object
            containing the `sample_metadata`, `feature_annotation` and the matrix for
            the given slice ranges.
        """
        _ssubset = self.get_sample_subset(sample_subset)
        _sample_indices = _ssubset.index.tolist()

        if not isinstance(feature_subset, (int, slice)):
            raise TypeError("feature indices must be continous; either a 'slice' or 'int' index.")
        _fsubset = self.get_feature_subset(feature_subset)
        start_findex = _fsubset["genarr_feature_start_index"].astype(int).min()
        end_findex = _fsubset["genarr_feature_end_index"].astype(int).max()

        # expand intervals
        final_rows = []
        for row in _fsubset.itertuples():
            for i, _ in enumerate(range(int(row.genarr_feature_start_index), int(row.genarr_feature_end_index))):
                final_rows.append(row._replace(starts=i + row.starts, ends=i + row.starts + 1))
        _feature_df = pd.DataFrame(final_rows)

        _msubset = self.get_matrix_subset((list(range(start_findex, end_findex)), _sample_indices))

        return GenomicArrayDatasetSlice(
            _ssubset,
            _feature_df,
            _msubset,
        )

    ####
    ## Dunder method to use `[]` operator.
    ####
    def __getitem__(
        self,
        args: Union[int, Sequence, tuple],
    ) -> GenomicArrayDatasetSlice:
        """Subset a ``GenomicArrayDataset``.

        Mostly an alias to :py:meth:`~.get_slice`.

        Args:
            args:
                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be extracted.

                Alternatively a tuple of length 1. The first entry specifies
                the rows (or cells) to retain based on their names or indices.

                Alternatively a tuple of length 2. The first entry specifies
                the rows (or cells) to retain, while the second entry specifies the
                columns (or features/genes) to retain, based on their names or indices.

        Raises:
            ValueError:
                If too many or too few slices provided.

        Returns:
            A :py:class:`~genomicarray.GenomicArrayDatasetSlice.GenomicArrayDatasetSlice` object
            containing the `sample_metadata`, `feature_annotation` and the matrix.
        """
        if isinstance(args, (str, int)):
            return self.get_slice(args, slice(None))

        if isinstance(args, tuple):
            if len(args) == 0:
                raise ValueError("At least one slicing argument must be provided.")

            if len(args) == 1:
                return self.get_slice(args[0], slice(None))
            elif len(args) == 2:
                return self.get_slice(args[0], args[1])
            else:
                raise ValueError(f"`{type(self).__name__}` only supports 2-dimensional slicing.")

        raise TypeError("args must be a sequence or a scalar integer or string or a tuple of atmost 2 values.")

    ####
    ## Misc methods.
    ####
    @property
    def shape(self):
        return (
            self._feature_annotation_tdb.nonempty_domain()[0][1] + 1,
            self._sample_metadata_tdb.nonempty_domain()[0][1] + 1,
        )

    def __len__(self):
        return self.shape[0]

    ####
    ## Printing.
    ####

    def __repr__(self) -> str:
        """
        Returns:
            A string representation.
        """
        output = f"{type(self).__name__}(number_of_rows={self.shape[0]}"
        output += f", number_of_columns={self.shape[1]}"
        output += ", at path=" + self._dataset_path

        output += ")"
        return output

    def __str__(self) -> str:
        """
        Returns:
            A pretty-printed string containing the contents of this object.
        """
        output = f"class: {type(self).__name__}\n"

        output += f"number_of_rows: {self.shape[0]}\n"
        output += f"number_of_columns: {self.shape[1]}\n"
        output += f"path: '{self._dataset_path}'\n"

        return output
