# MECH5961 — Assignment Section 8

Open this folder in MATLAB R2025a. Base MATLAB is sufficient.

```matlab
MECH5961_TMM_design  % Run the final design using the teacher's complete script.
run_submission       % Recompute required results, verification and all report plots.
make_report_figures  % Redraw report plots from the saved figure data only.
```

## Files to grade

| File | Assignment evidence |
|---|---|
| `MECH5961_TMM_design.m` | §8.5: final teacher script; only five student assignments changed |
| `data/design_parameters.csv` | Type, count, cavity depth, neck radius and extension in mm |
| `baseline/absorption.csv`, `final/absorption.csv` | §8.2/8.5: complete 200–2500 Hz curves, 1 Hz, 40 terms |
| `data/performance.csv` | Target-band mean, full-range peak frequency and peak absorption |
| `baseline/absorption.png`, `final/absorption.png` | Actual teacher Figure exports |
| `data/numerical_verification.csv` | §8.6: fixed final geometry at 1 Hz/40 and 0.25 Hz/80; signed mean change |
| `figures/` and `data/figure_data.mat` | Five main and three supplementary report plots, with their numeric source arrays |
| `results_handoff.md` | Writing facts, figure map, design rationale and model assumptions |

Both absorption CSVs have exactly `frequency_Hz,absorption`, with 2301 rows including 200 and 2500 Hz. Assignment §8.5 names `absorption.csv` but does not prescribe a header schema. The mean is trapezoidal integration over 800–1200 Hz divided by 400. The supplied teacher file has no CSV exporter; `run_submission` adds it externally. Each PNG is 600 dpi; `make_report_figures` reconstructs editable MATLAB plots from one MAT file.

## Search used to select the parameters

The original study screened 245 feasible cells on h=[10:2:20], r=[0.50:0.125:1.25], e=[0:2:10] mm, excluding e>=h. All 29,890 pairs used all 15 positive allocations. A diverse 24-cell subset supplied all 2,024 triples, each with 105 allocations; the subset then expanded to 36 cells. Bounded direct search used at least 8/10/12 initial starts for K=1/2/3, initial normalized step 0.10, unsuccessful-poll factor 0.5, h/e resolution 0.001 mm and r resolution 0.0001 mm. Counts were re-enumerated after every geometry change. Seeds were 20261021/20261022/20261023; promising runs received targeted continuation. `performance.csv` retains the selected single-, two- and three-type comparison. The MAT contains the plotted scan metrics and search traces, rather than full development logs.

The course common settings and nine original model functions are unchanged. The runner consolidates the reused model as local functions. The final design is best found in a nonconvex search; fixed-geometry count choices were exhaustive. Predictions assume independent sealed cavities under common surface pressure and a rigid, perfectly sealed backing. They do not establish physical performance or continuous global optimality.

The numeric study and supplied PNGs were computed in R2025a. CSV/MAT integrity and unchanged model source have been checked after consolidation. The consolidated runner has not been re-executed in the current restricted environment (MathWorks service error 5201).
