# Numerical and interchange contract

Updated local handout: `LAB2/MECH5961_Lab2_Acoustic_Source_Localisation (1).pdf`,
Sections 2.1–2.3 and 4.1. The obsolete 60 rev/s download is excluded. Course
experimental-methods slides support synchronized acquisition; sensor slides
explain sensitivity/phase limitations. Their proprietary content is not included.

MATLAB reference: 56 channels, fs=50000 Hz nominal, 500000 samples, 50 disjoint
10000-sample rectangular FFT blocks. Positive interior coefficients are
sqrt(2) FFT/L, DC/Nyquist unchanged. CSM is the mean of complex outer products;
all diagonal terms retained. With r the source-to-microphone distance,
T=exp(-i 2πfr/c0)/(4πr), w=T/sum(abs(T)^2), B=real(w' CSM w).
c0=343 m/s; h=1 m; no geometric calibration, channel reordering, smoothing,
diagonal removal or alternative beamformer. Numerical grid refinement changes
sampled scan positions, not the propagation model.

Every map is indexed [y,x,frequency,(macroblock)] in MATLAB and SciPy; no
transpose is applied. `geometry.mat` contains xy rows 1..56 in metres and
the monotonically increasing x/xb/xf axes (301/151/601 samples). Run MAT files
contain f (Hz), meanMaps (absolute nonnegative CBF power), macroMaps,
fineMaps/fineF, C and macroCSM (complex pressure cross spectra), and comparison
errors against prior validated 601 maps where available. These errors equal
zero for runs 1–3; run 0 has no prior 601 comparator and its zero initialized
field is not claimed as an independent baseline test. Python separately checks
scalar CBF values and grid coincidences for all runs.

`run3.mat` also stores spectrumFrequency (Hz) and spectrumPower (Pa² per bin).
Mic1 spectrum uses a periodic Hann window, energy normalization
FFT/sqrt(L sum(w²)), positive-bin RMS scaling and mean powers over 50 blocks.
Bin spacing 5 Hz; ENBW 7.5 Hz. This is a per-bin spectrum, not Pa²/Hz PSD.
SPL estimate=10log10(power/(20e-6 Pa)²). No detrending or additional sensitivity
coefficient is applied. The recorded Pa units do not certify calibration.

The prior independent spectral selection is rechecked: strongest locally
prominent tone in 1000–2000 and 2500–4000 Hz, using a 6 dB spectral contrast
to the local median (10<|Δf|≤50 Hz); in 4500–6500 Hz no such tone qualifies,
so the strongest measured bin is selected. This spectral selection contrast
is unrelated to spatial resolution. Frozen selected frequencies are 2000,
3000,6185 Hz. The last feature's mechanism is unassigned. Theory f_BPF=n fr
with n=2 and fr=100 rev/s gives 200 Hz and its second harmonic 400 Hz; guide
lines do not prove corresponding measured components.

Relative beamforming level is 10log10(B/max B) per entire 2D map, never Mic1
SPL. F11A samples y=0 without smoothing. F11D samples through the true dominant
grid peak; neither spatial axis is recentered. If quantified, lobe width is
the distance between the nearest strict flanking local minima, not half-power
width; missing minima give N/A. Nominal position error is Euclidean distance
to the confirmed nominal x/y, only for supported separate peaks.

Hashes and array dimensions are in frozen_manifest.csv and the numerical QA.
MATLAB v7 compression is lossless; SciPy reads original doubles. Frequency
arrays are explicitly cast to floating point before subtraction to prevent
unsigned integer wraparound from MATLAB's compact integer storage.
