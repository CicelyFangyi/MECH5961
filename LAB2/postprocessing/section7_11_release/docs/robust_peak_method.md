# Registered two-dimensional evidence rule

This release supersedes all prior fixed acoustic valley-depth resolution rules
and prior critical-frequency results. Neither half-power width nor an acoustic
dB cutoff determines the new status. The existing CBF propagation model is unchanged.

1. Detect all 8-connected maxima on nonnegative full 2D power before any source
   association. Normalize only to avoid amplitude-dependent floating tolerance.
   Ties are grouped at 64 double machine epsilons of the map maximum; reported
   power/levels retain the original values. Flat peak plateaus use the smallest
   row-major index deterministically. Flat shoulders with zero persistence are
   discarded. Any plateau touching the scan boundary is excluded as a source
   maximum and retained in the candidate audit.
2. Descending superlevel union-find constructs the elder-rule merge tree. A
   secondary maximum dies when its connected component joins a stronger one.
   A pair's saddle is its highest connecting bottleneck (common ancestor), not
   a minimum along the straight connecting line. ΔLs is the weaker peak level
   minus this saddle level, in dB; it is descriptive only. No valid pair means
   an empty/undefined value, never a fabricated zero.
3. Choose the strongest two interior maxima, then associate them left-to-right
   to independently confirmed nominal positions. Each must lie within half the
   source separation of its own nominal position. This preregistered geometric
   identity gate does not create, center, or split peaks. Other candidates,
   including boundary peaks, remain in T11_peak_candidates.csv.
4. Compare secondary-to-primary displacement to every single-source secondary
   maximum at the same frequency. A control conflict requires displacement
   agreement within max(40 mm, two grid cells) AND no greater saddle prominence
   than that control peak. This empirical conjunction flags sidelobe-like
   structure; it is not a universal acoustic threshold. Different source
   geometry limits the control's strength; control evidence alone is not proof
   of a spurious maximum. No candidate pair is accepted on geometry alone.
5. Use five synchronous consecutive 2 s whole-array macroblocks (10 rectangular
   non-overlapping 0.2 s FFT windows each). Their maps are calculated by MATLAB
   from intact 56-channel CSMs. Test all five blocks, all five leave-one-block-out
   means, and a 0–4 s versus 4–10 s split. Each must contain the same associated
   pair within 40 mm of the full-record estimates. Also report 30/50 mm identity
   sensitivity. This is a deterministic robustness screen, not a probability,
   confidence interval, or physical resolution distance. Five blocks from one
   record do not establish statistical independence; no bootstrap percentages
   or formal uncertainty intervals are reported.
6. Full-record maps use 301×301, block diagnostics use 151×151, both over
   [-0.5,0.5] m. The six prescribed conditions are independently computed on
   601×601; pair identity must agree. All stored maps use identical physics,
   coordinate orientation and normalization. The lower-resolution macroblock
   grid is explicitly an evidence limitation near shallow bifurcations.

## Frequency interpretation

The new scan is 1–6 kHz at 250 Hz steps, with 5 Hz FFT bins. A coarse sustained
candidate requires three consecutive sampled frequencies satisfying the
deterministic evidence rule. The first candidates are 1750 Hz (0.40 m) and
3750 Hz (0.20 m). These are **not published as f_onset**. In the 0.20 m record,
3250 Hz passes the 40 mm rule but not the 30 mm variant, and 3500 Hz loses block
support. Non-monotonic behavior prevents filling gaps by assumption.

Both onset results remain **indeterminate**: the coarse scan does not establish
the first sustained structure at finer frequency spacing, and the candidate
neighborhoods lack independent 601-grid refinement and block-grid convergence.
No numeric onset bracket is claimed from these incomplete neighborhoods. This
is a deliberately bounded result, not a claim of an intrinsic diffraction limit.
Further onset-specific analysis would require targeted finer-frequency and
block-grid checks; it is not a new physical experiment. Descriptive handout
scaling remains δ ∝ λ h / D, λ=c0/f, with no invented Rayleigh constant.

## Tests and limitations

Synthetic tests cover single/merged sources, separate peaks, plateau,
boundary maximum, known saddle, an off-line higher saddle path, a sidelobe
control conflict and missing block support. Synthetic values never enter
experimental tables. The conservative primary-two selection may reject a real
pair when a stronger unrelated maximum is present; complete candidate records
allow human review. No tuning was made to enforce monotonic resolution.
