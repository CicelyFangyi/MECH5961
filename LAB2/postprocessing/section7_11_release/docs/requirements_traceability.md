# V4 requirement traceability

| Requirement | Implementation/evidence | Artifact |
|---|---|---|
| Read-only course/model baseline | MATLAB rebuild, block_spectra, cbf_map; input hashes | input_manifest.csv; prior_matlab_baseline.txt |
| Actual geometry and unknown cases | geometry.mat; acquisition table; schematic-only panels a,d | F7A, F7B, T7_acquisition_geometry |
| Unconstrained 2D peaks and saddle persistence | peak_structure.py; candidate audit; synthetic tests | T11_peak_candidates; qa/peak_evidence.json |
| Sidelobe control | Run 0 same-frequency full maps; empirical displacement/prominence conjunction | qa/peak_evidence.json |
| Synchronous temporal evidence | MATLAB five intact-array macro-CSMs/maps; splits and leave-one-out | T11_frequency_evidence; F11B |
| Grid evidence | Independent 301/601 at prescribed frequencies; 151 block diagnostic grid declared | qa/acceptance.csv; T11_resolution_summary |
| No unsupported onset | Coarse consecutive-frequency scan and explicit missing finer-neighborhood evidence | T11_onset_summary: indeterminate |
| Horizontal profiles and supported coordinates | Fine maps, y=0, projected peak marks, full xy in tables | F11A; T11_resolution_summary |
| Mic1 spectral identity and BPF | Hann energy-normalized single-channel power, theoretical 200/400 Hz guides | F11C; T11_propeller_frequency_summary |
| Actual x/y peak profiles and data-defined spans | Dominant 2D peak cuts, strict flanking local minima | F11D; T11_propeller_frequency_summary |
| Conditional nominal accuracy | Errors only for supported pair; optional figure | F11E; T11_position_accuracy |
| Unified Python art and dimensions | render_release.py; contracts; physical PDF/PNG/text/alignment checks | Seven PNG/SVG/PDF triples; qa |
| Frozen independent reproduction | Five derived MAT files; clean venv; no private inputs | render_release.py --from-frozen; clean_reproduction.json |
| Optional raw rebuild | Project-owned MATLAB only; external authorized inputs | code/matlab/rebuild_release_data.m |
| Single authoritative output folder | All figures and tables together | figures_tables/FIGURE_TABLE_MANIFEST.csv |
| Scope and publication | No Section 8 files, source raw data or teacher code; isolated branch | github_release_receipt.md |

Indeterminate onset is an allowed undefined scientific result, not a hidden
substitution of prior critical frequencies. Further onset-specific computation
is explicitly outside the evidence claimed by this release.
