# Section 7 / Appendix 11 release — V4

Status: HUMAN_REVIEW_REQUIRED. No report prose and no Section 8 artwork included.
Authoritative figures and tables are exclusively in `../figures_tables/`.
Seven figures are provided individually as 600 dpi PNG, editable-text SVG and
embedded-font PDF. F11E is an optional review candidate, not automatic report inclusion.

## Reproduce without private recordings or MATLAB

From this release directory (Python 3.13 tested):

```text
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python code/python/test_peak_structure.py
.venv\Scripts\python code/python/render_release.py --from-frozen
.venv\Scripts\python code/python/validate_release.py
```

On Linux/macOS use `.venv/bin/python`. Entry scripts resolve the release root
from their own location; any current working directory is supported. No skill
installation is needed. Arial is preferred when installed, with DejaVu Sans
fallback; fonts are not redistributed. Renderer/font metadata can differ across
platforms; numerical tables and physical layouts remain the contract.

## Optional authorized raw-data rebuild

In MATLAB add `code/matlab` to the path, then call
`rebuild_release_data('/path/to/local/LAB2')`. The input directory must contain
`Beamforming/Coor_56Array.txt` and the four converted `MicData_0..3.mat` in
`Group2-20261116T100220Z-1-001`. The original TDMS conversion is intentionally
not redistributed. Existing frozen files are retained; move them outside the
release explicitly when requesting a full rebuild. A public clone cannot
rebuild raw measurements without these authorized external inputs.

The included MATLAB helpers are project-owned implementations of the audited
CBF equations. No teacher script, handout, original coordinate text, raw audio
or raw measurement archive is included. MATLAB writes numerical MAT arrays
only; Python creates all release artwork. The plotting entry recalculates all
tables and peak evidence from those arrays rather than embedding figure values.

See `FIGURE_CONTRACTS.md`, `robust_peak_method.md`, `methods_and_definitions.md`,
`provenance_and_limits.md`, and `validation_report.md`. Full figure titles,
external-caption limitations, placement and source links are in
`figures_tables/FIGURE_TABLE_MANIFEST.csv`.
