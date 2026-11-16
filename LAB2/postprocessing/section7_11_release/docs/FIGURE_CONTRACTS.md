# Figure contracts (registered before rendering)

Python/Matplotlib exclusively; unmodified MATLAB numerical maps. All panels use
8–9 pt essential text, 9 pt axis labels, 1.4 pt data lines, 0.8 pt spines,
white backgrounds, blue/ochre/charcoal. Labels below panels override the skill's
generic upper-left label suggestion. PNG 600 dpi, editable SVG, embedded-font PDF.
Exact physical export size is retained (no tight-bbox page resizing).

| Figure | Scientific purpose / evidence | Architecture and width × height | Annotation risk |
|---|---|---|---|
| F7A | Preserve the measured array layout and coordinate convention | Single quantitative panel, 100 × 95 mm | Equal aspect; no origin text or inferred microphone IDs |
| F7B | Distinguish four experiments and actual asymmetric speaker arrangement | 2×2 setup schematic, 183 × 145 mm | Single/propeller qualitative panels have no numeric coordinates; quantitative panels share axes |
| F11A | Compare horizontal profiles across spacing and frequency without manufacturing peaks | 2×3 quantitative grid, 183 × 165 mm | Only supported 2D peaks projected to y=0; saddle descriptor is not a linecut valley |
| F11B | Display deterministic 2D peak evidence across analysis frequencies | Two equal panels, 183 × 90 mm | No support percentages, fixed acoustic cutoffs or unsupported onset |
| F11C | Distinguish measured Mic1 spectrum from theoretical BPF guides | Two complementary panels, 183 × 90 mm | Calibration caveat external; selected 6185 Hz feature has no assigned mechanism |
| F11D | Describe actual propeller map profiles in both directions | 2×3 quantitative grid, 183 × 165 mm | No recentering or half-power width; local-minimum spans may be undefined |
| F11E | Assess nominal position discrepancy only for supported pairs | Single quantitative panel, 120 × 95 mm | Review candidate; physical wiring/geometry uncertainty prevents certified accuracy |

No figure supports repeated-experiment confidence intervals. Time blocks come
from one 10 s recording. F11B shows deterministic evidence stages, not an
acoustic resolution scale. External captions and full numerical values reside
in the figure/table manifest and tables, not paragraphs inside the artwork.
