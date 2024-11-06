"""A dataloader using TileDB files in the pytorch-lightning framework.

This class provides a dataloader using the generated TileDB files built using the
:py:func:`~genomicarrays.build_genomicarray.build_genomicarray`.

Example:

    .. code-block:: python

        from genomicarrays.dataloader import TorchDataset

        ds = TorchDataset(
            dataset_path="/path/to/genarr/dir"
        )

        print(ds[0])
"""

import logging
from typing import Union

from torch.utils.data import Dataset

from .GenomicArrayDataset import GenomicArrayDataset

log = logging.getLogger(__name__)


__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


class TorchDataset(Dataset):
    """A class that extends pytorch :py:class:`~torch.utils.data.Dataset` to enumerate features and samples using
    TileDB."""

    def __init__(self, dataset_path: Union[str, GenomicArrayDataset]):
        """Initialize a ``gaTorchDataset``.

        Args:
            dataset_path:
                Path to the directory containing the TileDB files.

                Alternatively, may also provide a ``GenomicArrayDataset``
                object.
        """
        self._dataset = dataset_path
        if isinstance(dataset_path, str):
            self._dataset = GenomicArrayDataset(dataset_path=dataset_path)

    def __len__(self):
        return len(self._dataset)

    def __getitem__(self, idx):
        return self._dataset[idx, :]

    def __repr__(self) -> str:
        """
        Returns:
            A string representation.
        """
        output = f"{type(self).__name__}("
        output += f"number_of_features={self._dataset.shape[0]}"
        output += f"number_of_samples={self._dataset.shape[1]}"
        output += ")"

        return output

    def __str__(self) -> str:
        """
        Returns:
            A pretty-printed string containing the contents of this object.
        """
        output = f"class: {type(self).__name__}\n"
        output += f"number_of_features: {self._dataset.shape[0]}\n"
        output += f"number_of_samples: {self._dataset.shape[1]}\n"

        return output
