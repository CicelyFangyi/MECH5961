# MECH5961 Lab 2 — Acoustic Source Localisation

## 1. Prerequisites

Python 3.13 (64-bit), approximately 2 GB available RAM, and locally installed
Times New Roman regular and italic fonts. Missing fonts cause an explicit failure;
no font substitution is permitted. Install dependencies once:

```text
python -m pip install -r postprocessing/requirements.txt
```

The program does not require MATLAB, a plotting skill, or saved development maps.
Runtime depends on CPU speed; reconstruction evaluates the original grids and
all required synchronized time-block maps, rather than loading cached results.

## 2. Input data and folder structure

```text
LAB2/
  Beamforming/Coor_56Array.txt
  Group2-20261116T100220Z-1-001/MicData_0.mat ... MicData_3.mat
  postprocessing/run_lab2.py
  postprocessing/requirements.txt
  postprocessing/README.md
  figures_tables/
```

The local TA package retains the complete original input directories, including
teacher MATLAB programs and original TDMS files. The program reads only the
comma-delimited 56-by-2 coordinate file and the four converted MAT records.
Each record must contain `ConvertedData.Data.MeasuredData`, with 56 ordered
`Untitled/Sound Pressure_0` through `Sound Pressure_55` channels, 500000 finite
samples per channel, Pascals metadata, and approximately 50000 Hz acquisition.

**Public GitHub distribution:** raw MAT/TDMS recordings and teacher files are
not included. Obtain these authorized course/Group 2 inputs from the teaching
team or experiment group and place them in the folders above. A public clone
alone can inspect the delivered figures and code but cannot reconstruct the
measurement results without these external inputs. No LFS availability or
permission to publicly distribute course materials is assumed.

## 3. Run

```text
python postprocessing/run_lab2.py
```

An absolute script path works from any current directory. LAB2 is resolved from
the script location. An optional `--output PATH` selects another output folder.
The default is `LAB2/figures_tables`. Inputs are processed one record at a time.
All intermediate maps exist only within the current run. Outputs are prepared
in a temporary directory and replace the final output directory only after
checks pass. Unknown files in an existing output folder block replacement.
Failures report the exception and preserve the previous completed output.

## 4. Outputs

Seven figures, each as 600 dpi PNG, editable SVG and embedded-font vector PDF:

- F7A: actual 56-channel array geometry.
- F7B: four experimental configurations.
- F11A: six two-source horizontal profiles.
- F11B: two-dimensional saddle evidence and temporal stability.
- F11C: microphone-1 propeller spectrum and BPF detail.
- F11D: peak-through propeller profiles and flanking-minimum spans.
- F11E: nominal positions and signed localisation-error components.

CSV tables describe acquisition/geometry, six-condition resolution, frequency
evidence, onset limitations, propeller features, localisation errors and profile
boundaries. `FIGURE_TABLE_MANIFEST.csv` supplies titles, units, methods, limitations,
report placements and figure sizes. Empty numeric fields mean not defined, not
zero. No Section 8 figures, candidate-peak dumps or report prose are generated.

## 5. Scientific method

Run 0 is a single loudspeaker; run 1 has nominal x positions (-0.20, 0.00) m;
run 2 has (-0.20, +0.20) m; run 3 is a two-bladed propeller. Array coordinates
are in the z=0 plane; nominal sources lie in z=1 m, with y=0 for the dual sources.
F7B draws the single source and propeller on the z axis as requested schematics,
not independently measured absolute positions.

CBF uses 50 non-overlapping 10000-sample rectangular blocks at 50000 Hz, with
5 Hz FFT bins and one-sided RMS amplitudes `sqrt(2)*FFT/N` at the positive
analysis frequencies. CSM is the average outer product across blocks. The full
CSM, including its diagonal, is retained. For range r, the Green function is
`T=exp(-i*2*pi*f*r/343)/(4*pi*r)` and unit-response weights are
`w=T/sum(abs(T)^2)`. Beamforming power is `real(w.conj().T @ CSM @ w)`.
Map levels use `10*log10(B/max(B))` over the full map, not microphone SPL.

The scan is x,y in [-0.5,+0.5] m. Mean maps use 301 points per dimension;
five synchronized 2-second macroblocks use 151; final 1/3/6 kHz maps and the
three propeller maps use 601. The frequency evidence scan is 1000:250:6000 Hz.
The propeller microphone spectrum uses periodic Hann windows with window-energy
normalization, one-sided bin power and 50 block averages. Selected frequencies
are 2000, 3000 and 6185 Hz, reconfirmed by the recorded spectral-selection rule.
Theoretical BPF=2 blades * 100 revolutions/s=200 Hz; 400 Hz is its second harmonic.
These theoretical markers are distinguished from measured local spectral peaks.

Two-dimensional superlevel connectivity supplies peak/saddle prominence.
Detected peaks are associated with nominal source coordinates only after map
peak detection. Single-source sidelobes, five synchronized blocks, leave-one-out
and split-block consistency, 30/40/50 mm identity-tolerance sensitivity and
301/601-grid agreement provide evidence. There is no fixed 3 dB resolution
threshold. Regression states are: 0.40 m unresolved at 1 kHz and robust at
3/6 kHz; 0.20 m unresolved at 1/3 kHz and robust at 6 kHz. Onset frequency
remains indeterminate. Coarse candidates are diagnostic only. Main-lobe spans
use the nearest strict minimum on each side; missing boundaries remain undefined.
Errors are ex=xhat-xnom, ey=yhat-ynom and e=hypot(ex,ey).

## 6. Known measurement limitations

Run mapping and dual-source coordinates derive from participant recollection.
Pascals metadata alone does not establish the microphone sensitivity-calibration
chain; SPL is therefore an estimate referenced to 20 microPa. Physical channel
wiring to coordinate-file rows is unverified. Nominal position residuals are not
certified metrological accuracy. Time blocks are not independent experiments;
no unsupported statistical error bars or physical critical-frequency claims
are supplied. Source glyph sizes are schematic. F11E remains subject to human
review for inclusion. The physical CBF model and original data are unchanged.
