# Verified GitHub release receipt

Status: **HUMAN_REVIEW_REQUIRED**. Date: 2026-11-16.

- Repository: https://github.com/CicelyFangyi/MECH5961
- Release branch: https://github.com/CicelyFangyi/MECH5961/tree/lab2/robust-figures-20261116/LAB2/postprocessing/section7_11_release
- Base main SHA, verified unchanged after payload push:
  `3e5570d366892c8e18c174ff38f0251f0aa5f64c`
- Scientific payload commit:
  `1d3009a1799b440643ed267c4f3dd579df22cd62`
- This receipt and remote-verification records are added in a follow-up commit.
  The final branch-tip SHA is verified after that push and included in the
  delivery message (a receipt cannot contain its own Git commit hash).

GitHub's recursive remote tree returned all 101 payload files, without
truncation. Each remote Git blob SHA was independently compared with the
local file bytes; all matched. This includes every PNG/SVG/PDF/CSV, Python
entry, MATLAB helper and all five frozen MAT files. See
`qa/remote_payload_tree.json` for the server-returned tree.

A separate fresh clone of the published branch, using the clean virtual
environment, ran `render_release.py --from-frozen` and `validate_release.py`.
No MATLAB or private raw inputs were used. All eight CSVs are byte-identical
and all seven PNGs pixel-identical to the local release. All 23 numerical and
export checks passed, with maximum relative discrepancy 1.128574536526135e-15.
See `qa/remote_clone_reproduction.json`.

Representative verified payload links:

- [F11A PNG](https://github.com/CicelyFangyi/MECH5961/blob/1d3009a1799b440643ed267c4f3dd579df22cd62/LAB2/postprocessing/section7_11_release/figures_tables/F11A_two_source_profiles.png)
- [F11A SVG](https://github.com/CicelyFangyi/MECH5961/blob/1d3009a1799b440643ed267c4f3dd579df22cd62/LAB2/postprocessing/section7_11_release/figures_tables/F11A_two_source_profiles.svg)
- [F11A PDF](https://github.com/CicelyFangyi/MECH5961/blob/1d3009a1799b440643ed267c4f3dd579df22cd62/LAB2/postprocessing/section7_11_release/figures_tables/F11A_two_source_profiles.pdf)
- [Resolution table](https://github.com/CicelyFangyi/MECH5961/blob/1d3009a1799b440643ed267c4f3dd579df22cd62/LAB2/postprocessing/section7_11_release/figures_tables/T11_resolution_summary.csv)
- [Python entry](https://github.com/CicelyFangyi/MECH5961/blob/1d3009a1799b440643ed267c4f3dd579df22cd62/LAB2/postprocessing/section7_11_release/code/python/render_release.py)
- [MATLAB rebuild entry](https://github.com/CicelyFangyi/MECH5961/blob/1d3009a1799b440643ed267c4f3dd579df22cd62/LAB2/postprocessing/section7_11_release/code/matlab/rebuild_release_data.m)

Only `LAB2/postprocessing/section7_11_release/` differs from main. Original
teacher sources, recordings, acquired audio, temporary caches, virtual
environment and Section 8 figures are absent from the added tree. No force
push, PR, merge, or modification of main/Lab 1 was performed. The publication
checkout is clean after the receipt commit; regenerating SVGs in the separate
verification clone can change harmless renderer IDs/metadata.
