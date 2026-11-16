# V4 acceptance record

Stop state: **HUMAN_REVIEW_REQUIRED**. This is a methods/QA record, not report prose.

## Gate A — source and MATLAB reference: PASS

All 47 entries in the prior source inventory were rehashed unchanged. The
updated local handout and actual supplied CBF implementation were reread.
Teacher files and original measurements were not modified or redistributed.
The prior independently written FFT/CBF baseline reports zero CSM discrepancy
and 1.6823e-16 relative map discrepancy. This run's optional raw-data rebuild
read the converted MAT recordings once per run, preserving all 56 channels,
and exported synchronous macroblock maps without TDMS reconversion.
Runs 1–3 prescribed 601-grid maps match the previous validated cache exactly.

The 153.2 MB frozen numerical bundle retains double-precision full maps,
macroblock maps and compact CSMs; it is derived data, not acquired audio. Each
file is below 50 MB. A field-level shape/dtype/unit schema, SHA-256 manifest
and input hashes are supplied. No original MicData MAT/TDMS is needed to render.

## Gate B — bounded scientific result: PASS with onset indeterminate

Six synthetic test methods cover eight failure/edge cases, including an
off-line connecting saddle and absent time-block support. All passed.
Unconstrained maxima and 0D merge persistence precede source association.
All peaks, boundary flags, control evidence, five macroblocks, five leave-one-
out tests, two temporal splits and identity-tolerance sensitivity are retained.
No fixed acoustic dB threshold, bootstrap percentage or invented interval is used.

| Separation | 1 kHz | 3 kHz | 6 kHz |
|---|---|---|---|
| 0.40 m | UNRESOLVED | ROBUST_DUAL_PEAK | ROBUST_DUAL_PEAK |
| 0.20 m | UNRESOLVED | UNRESOLVED | ROBUST_DUAL_PEAK |

Supported 601-grid 2D saddle depths are respectively 15.6609451577,
17.5655727885 and 11.0311158861 dB; these describe structure and do not define
the status. At 0.20 m / 3 kHz only one of five macroblocks supports the same
pair, so no two-source positions/errors are reported. At 3250 Hz the result
depends on identity tolerance; at 3500 Hz block support is lost. The coarse
candidate values 1750/3750 Hz are explicitly not f_onset. Both onset cells are
indeterminate because finer frequency and candidate-neighborhood grid evidence
are missing. This limit is exposed rather than relabeling old threshold crossings.

## Gate C — numerical and figure integrity: PASS

23 interchange/export checks passed. Maximum independent scalar/grid relative
discrepancy: 1.128574536526135e-15. Synchronous macroblock means agree with full
maps at coincident pixels; 301/601 common pixels and four scalar steering
evaluations per fine map agree to floating precision. The original measured
spectrum selects 2000,3000,6185 Hz again. Tables are recalculated from these
arrays, including peak positions, descriptive saddle depths and local-minimum
spans. Boundary-truncated spans are empty rather than extrapolated.

All seven PNG/SVG/PDF triples reopen. PNG dimensions match 600 dpi and PDF
page sizes match the declared 100/120/183 mm widths. SVG retains text and PDF
embeds fonts. Axis labels are 9 pt, ticks 8 pt, essential annotations 8 pt;
small mathematical subscripts reach 5.6 pt, above the 5 pt glyph floor.
All comparable panel plot areas pass the 1.5 pt alignment gate. Every final
PDF has zero collision FAILs and zero WARNs. F11A mathematical accent/fill
overlays are intentional, legible text glyphs rather than curve masking.

Individual final-size previews and every panel were visually inspected:
axes, labels below panels, peak projections, neutral reference guides, units,
color/marker distinctions and legends are readable without upper titles.
F7A has no origin coordinate text; F7B contains all four experiments, with
qualitative a,d separated from quantitative b,c. F11E remains optional.
There are no Section 8 artwork assets in the release.

Source-preflight warnings are reviewed: TIFF is not required by V4 (600 dpi
PNG is supplied); dynamic figure widths are tested in the actual PDF rather
than inferred statically; logarithms have positive-data assertions or explicit
machine-tiny guards. No warning changes scientific data.

## Gate D — clean standalone reproduction: PASS

A new Python 3.13 virtual environment installed the pinned requirements.
The release was copied into an isolated checkout containing no private raw
inputs. The documented entry and tests succeeded without MATLAB, TDMS or
the local skill. All CSV bytes and PNG pixels match across the original and
clean reproduction. SVG internal IDs/timestamps and renderer metadata may
differ; the physical layouts and data are identical. See clean_reproduction.json,
environment_lock.txt and acceptance.csv in qa. Full MATLAB rebuild is optional
and explicitly requires separately authorized local recordings.

## Gate E — publication

See github_release_receipt.md for the verified branch, commit and remote tree.
Only the complete release directory is added at its identical repository path;
main and Lab 1 are not edited. No PR or merge is requested or created.

Remaining human review: whether to include optional F11E; acceptance of nominal
geometry and unverified physical wiring; calibration limits on absolute SPL;
and whether further onset-specific analysis is needed before report writing.
