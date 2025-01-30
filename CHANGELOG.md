# Changelog

## Version 0.2.1 - 0.2.2

- The aggregate function is expected to return either a scalar value or a 1-dimensional NumPy ndarray. If the later, users need to specify the expected dimension of the summarization. All values will be flattenned eventually.
- Remove expanding the intervals to conform to output length; this is now incompatible with coercions to anndata and summarized experiments and has been removed.

## Version 0.2.0

- chore: Remove Python 3.8 (EOL)
- precommit: Replace docformatter with ruff's formatter

## Version 0.1.0

- Initial version of the package
