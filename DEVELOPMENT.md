# mafutils Development Notes

- `tests/example.maf.scaffold.idx` is preserved as an older scaffold-index
  fixture for comparison.
- `tests/example.maf.scaffold.regenerated.idx` is the scaffold index produced
  by the current `mafutils index` implementation.
- The real production scaffold indexes checked so far match the regenerated
  convention, not the preserved older fixture.

## Releasing to PyPI

- Publishing is triggered by pushing a Git tag that matches `v*` (for example
  `v0.1.1`).
- Ensure the GitHub repository secret `PYPI_API_TOKEN` is set before releasing.
- Commit release-related changes first, then create and push the tag:

```bash
git add .github/workflows/publish-pypi.yml pyproject.toml README.md DEVELOPMENT.md
git commit -m "release: prepare v0.1.1"
git tag -a v0.1.1 -m "v0.1.1"
git push origin main v0.1.1
```
