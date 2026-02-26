# DALES case generator docs

This folder contains Sphinx configuration for generating HTML documentation
for the DALES case generator, based on Python docstrings.

## Building the docs

From the project root (where `docs/` lives) and with `sphinx` installed in the
`phd` environment, run:

```bash
cd docs
sphinx-build -b html . _build/html
```

Then open `_build/html/index.html` in a browser.

If you prefer `make`:

```bash
cd docs
make html
```

## Extending the docs

- Add more modules or classes to `docs/api.rst` under the `.. autosummary::` list.
- Write narrative pages (e.g. `usage.rst`, `cases.rst`) and include them in the
  `.. toctree::` in `docs/index.rst`.
- Keep docstrings in the codebase up to date; Sphinx will pull those into the
  rendered HTML automatically.
