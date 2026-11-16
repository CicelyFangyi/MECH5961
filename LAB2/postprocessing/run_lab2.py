"""MECH5961 Lab 2 Acoustic Source Localisation.

Read Group 2 ConvertedData MAT pressure records and the supplied array coordinates.
Run: python postprocessing/run_lab2.py (from any working directory).
Output: seven 600 dpi PNG/SVG/PDF figures and formal CSV tables in LAB2/figures_tables.
The teacher conventional-beamforming physical model is preserved. Absolute
pressure calibration and physical channel wiring remain unverified. No raw
input or teacher source is modified. Runtime intermediates are temporary only.
"""

from pathlib import Path
import csv, json, os, sys, tempfile, shutil, hashlib, argparse

os.environ["OPENBLAS_NUM_THREADS"] = "2"
os.environ["OMP_NUM_THREADS"] = "2"
try:
    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.lines import Line2D
    from scipy.io import loadmat
    from scipy import ndimage
    from scipy.spatial import ConvexHull
    from numba import njit
    import pymupdf
except ImportError as error:
    raise SystemExit(
        "Missing dependency: "
        + str(error)
        + ". Install postprocessing/requirements.txt."
    )


# 1. Raw-input reconstruction and teacher CBF
def cbf_power(csm, xy, axis, frequency):
    """Teacher CBF: T=exp(-ikr)/(4*pi*r), w=T/(T* T), B=Re(w* C w).

    Full CSM, including its diagonal; no smoothing or channel permutation.
    Chunking changes memory use only. Coordinates and propagation distance: m.
    """
    xx, yy = np.meshgrid(axis, axis)
    output = np.empty(xx.size)
    for first in range(0, xx.size, 8192):
        selection = slice(first, first + 8192)
        distance = np.sqrt(
            (xy[:, 0, None] - xx.ravel()[selection]) ** 2
            + (xy[:, 1, None] - yy.ravel()[selection]) ** 2
            + 1.0
        )
        transfer = np.exp(-2j * np.pi * frequency * distance / 343) / (
            4 * np.pi * distance
        )
        weights = transfer / np.sum(np.abs(transfer) ** 2, axis=0)
        output[selection] = np.real(np.sum(weights.conj() * (csm @ weights), axis=0))
    assert output.min() >= -1e-10 * max(1, np.abs(output).max())
    return np.maximum(output, 0).reshape(xx.shape)


def reconstruct(lab_root):
    """Read one complete MAT record at a time; retain only spectra and maps."""
    xy = np.loadtxt(lab_root / "Beamforming" / "Coor_56Array.txt", delimiter=",")
    assert xy.shape == (56, 2) and np.isfinite(xy).all()
    geometry = dict(
        xy=xy,
        x=np.linspace(-0.5, 0.5, 301),
        xb=np.linspace(-0.5, 0.5, 151),
        xf=np.linspace(-0.5, 0.5, 601),
        frequencies=np.arange(1000, 6001, 250),
        fs=50000,
        L=10000,
        c0=343,
        h=1,
    )
    records = []
    for run in range(4):
        path = lab_root / "Group2-20261116T100220Z-1-001" / f"MicData_{run}.mat"
        print(f"Reading {path.name}", flush=True)
        channels = loadmat(path, simplify_cells=True)["ConvertedData"]["Data"][
            "MeasuredData"
        ]
        frequencies = (
            geometry["frequencies"] if run < 3 else np.array([2000, 3000, 6185])
        )
        spectra = np.empty((56, 50, len(frequencies)), complex)
        metadata_rates = []
        for channel in range(56):
            item = channels[channel + 2]
            assert item["Name"] == f"Untitled/Sound Pressure_{channel}", (path, channel)
            properties = {prop["Name"]: prop["Value"] for prop in item["Property"]}
            rate = 1 / float(properties["wf_increment"])
            assert abs(rate - 50000) < 0.001 and properties["unit_string"] == "Pascals"
            metadata_rates.append(rate)
            pressure = np.asarray(item["Data"], dtype=np.float64)
            assert pressure.shape == (500000,) and np.isfinite(pressure).all()
            # MATLAB reshape(L,50) is reproduced by contiguous Python blocks.
            blocks = pressure.reshape(50, 10000)
            transformed = np.fft.rfft(blocks, axis=1) / 10000
            spectra[channel] = np.sqrt(2) * transformed[:, frequencies // 5]
            if run == 3 and channel == 0:
                window = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(10000) / 10000)
                windowed = np.fft.rfft(blocks * window, axis=1) / np.sqrt(
                    10000 * np.sum(window**2)
                )
                power = np.mean(np.abs(windowed) ** 2, axis=0)
                power[1:-1] *= 2  # one-sided bin power, DC/Nyquist undoubled
        del channels, item, pressure, blocks, transformed
        fine_frequencies = (
            frequencies[np.isin(frequencies, [1000, 3000, 6000])]
            if run < 3
            else frequencies
        )
        record = dict(
            run=run,
            f=frequencies,
            fineF=fine_frequencies,
            metadata_fs_Hz=float(np.mean(metadata_rates)),
            C=np.empty((56, 56, len(frequencies)), complex),
            macroCSM=np.empty((56, 56, len(frequencies), 5), complex),
            meanMaps=np.empty((301, 301, len(frequencies))),
            macroMaps=np.empty((151, 151, len(frequencies), 5)),
            fineMaps=np.empty((601, 601, len(fine_frequencies))),
        )
        if run == 3:
            record.update(spectrumFrequency=np.arange(5001) * 5, spectrumPower=power)
        for index, frequency in enumerate(frequencies):
            amplitudes = spectra[:, :, index]
            csm = amplitudes @ amplitudes.conj().T / 50
            record["C"][:, :, index] = csm
            record["meanMaps"][:, :, index] = cbf_power(
                csm, xy, geometry["x"], frequency
            )
            for block in range(5):
                subset = amplitudes[:, block * 10 : (block + 1) * 10]
                macro_csm = subset @ subset.conj().T / 10
                record["macroCSM"][:, :, index, block] = macro_csm
                record["macroMaps"][:, :, index, block] = cbf_power(
                    macro_csm, xy, geometry["xb"], frequency
                )
            if frequency in fine_frequencies:
                fine_index = int(np.flatnonzero(fine_frequencies == frequency)[0])
                record["fineMaps"][:, :, fine_index] = cbf_power(
                    csm, xy, geometry["xf"], frequency
                )
            print(f"Run {run}: {frequency} Hz complete", flush=True)
        records.append(record)
    return geometry, records


# 2. Two-dimensional peak structure
TIE = 64 * np.finfo(float).eps
LOCATION_TOL_M = 0.04  # identity tolerance, not a resolution distance


@njit(cache=False)
def merge_tree(v, ny, nx):
    order = np.argsort(-v, kind="mergesort")
    uf = np.full(v.size, -1, np.int64)
    elder = np.arange(v.size)
    parent = np.full(v.size, -1, np.int64)
    death = np.full(v.size, np.nan)
    for p in order:
        uf[p] = p
        iy, ix = p // nx, p % nx
        for dy in range(-1, 2):
            for dx in range(-1, 2):
                yy, xx = iy + dy, ix + dx
                if (dy == 0 and dx == 0) or yy < 0 or yy >= ny or xx < 0 or xx >= nx:
                    continue
                q = yy * nx + xx
                if uf[q] < 0:
                    continue
                a, b = p, q
                while uf[a] != a:
                    uf[a] = uf[uf[a]]
                    a = uf[a]
                while uf[b] != b:
                    uf[b] = uf[uf[b]]
                    b = uf[b]
                if a == b:
                    continue
                ea, eb = elder[a], elder[b]
                if v[ea] < v[eb] or (v[ea] == v[eb] and ea > eb):
                    a, b, ea, eb = b, a, eb, ea
                parent[eb] = ea
                death[eb] = v[p]
                uf[b] = a
                elder[a] = ea
    return parent, death


def analyse(power, axis, reference=None):
    p = np.asarray(power, float)
    assert p.ndim == 2 and p.shape == (len(axis), len(axis))
    assert np.all(np.isfinite(p)) and p.min() >= 0 and p.max() > 0
    # Quantization only groups floating ties for topology; reported values use raw B.
    v = np.rint(p / p.max() / TIE) * TIE
    candidate = v == ndimage.maximum_filter(v, size=3, mode="nearest")
    labels, n = ndimage.label(candidate, np.ones((3, 3)))
    parent, death = merge_tree(v.ravel(), *v.shape)
    peaks = []
    for k in range(1, n + 1):
        yy, xx = np.where(labels == k)
        index = int(np.min(yy * len(axis) + xx))
        y, x = divmod(index, len(axis))
        boundary = bool(
            np.any(
                (xx == 0) | (xx == len(axis) - 1) | (yy == 0) | (yy == len(axis) - 1)
            )
        )
        saddle = death[index]
        if np.isfinite(saddle) and v[y, x] - saddle <= TIE:
            continue  # Flat shoulders are not separate maxima.
        peaks.append(
            dict(
                index=index,
                x_m=float(axis[x]),
                y_m=float(axis[y]),
                level_dB=float(10 * np.log10(p[y, x] / p.max())),
                prominence_dB=(
                    float(10 * np.log10(v[y, x] / saddle)) if saddle > 0 else None
                ),
                boundary=boundary,
                plateau_cells=len(xx),
            )
        )
    peaks.sort(key=lambda z: (-z["level_dB"], z["index"]))
    interior = [z for z in peaks if not z["boundary"]]
    pair = None
    depth = None
    saddle = None
    if reference is not None and len(interior) >= 2:
        # Primary pair is the strongest two interior maxima, chosen BEFORE matching.
        a, b = sorted(interior[:2], key=lambda z: z["x_m"])
        ref = np.array([[reference[0], 0], [reference[1], 0]])
        coords = np.array([[a["x_m"], a["y_m"]], [b["x_m"], b["y_m"]]])
        if np.all(
            np.linalg.norm(coords - ref, axis=1) < (reference[1] - reference[0]) / 2
        ):
            pair = [a, b]

            def chain(index):
                d = {index: np.inf}
                bottleneck = np.inf
                while parent[index] >= 0:
                    bottleneck = min(bottleneck, death[index])
                    index = int(parent[index])
                    d[index] = bottleneck
                return d

            ca, cb = chain(a["index"]), chain(b["index"])
            saddle = max(min(ca[z], cb[z]) for z in ca.keys() & cb.keys())
            depth = float(min(a["level_dB"], b["level_dB"]) - 10 * np.log10(saddle))
    return dict(
        peaks=peaks,
        pair=pair,
        saddle_depth_dB=depth,
        saddle_relative_power=None if saddle is None else float(saddle),
    )


def same_pair(a, b, tol=LOCATION_TOL_M):
    if a["pair"] is None or b["pair"] is None:
        return False
    return all(
        np.hypot(x["x_m"] - y["x_m"], x["y_m"] - y["y_m"]) <= tol
        for x, y in zip(a["pair"], b["pair"])
    )


def sidelobe_conflict(result, control, grid_step):
    """Control flags matching secondary displacement AND no stronger persistence.
    Different control geometry prevents treating this alone as proof of an artifact.
    """
    if result["pair"] is None:
        return False
    p = sorted(result["pair"], key=lambda z: -z["level_dB"])
    c = [z for z in control["peaks"] if not z["boundary"]]
    if len(c) < 2:
        return False
    delta = np.array([p[1]["x_m"] - p[0]["x_m"], p[1]["y_m"] - p[0]["y_m"]])
    for q in c[1:]:
        d = np.array([q["x_m"] - c[0]["x_m"], q["y_m"] - c[0]["y_m"]])
        if np.linalg.norm(delta - d) <= max(LOCATION_TOL_M, 2 * grid_step):
            if (
                q["prominence_dB"] is not None
                and result["saddle_depth_dB"] <= q["prominence_dB"]
            ):
                return True
    return False


def minima_span(axis, profile, index):
    ids = (
        np.flatnonzero((profile[1:-1] < profile[:-2]) & (profile[1:-1] < profile[2:]))
        + 1
    )
    left = ids[ids < index]
    right = ids[ids > index]
    if not len(left) or not len(right):
        return None
    return float(axis[right[0]] - axis[left[-1]])


# 3. Evidence and physical metrics
def csvwrite(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def process(root, G, runs):
    out = root / "figures_tables"
    qa = root / "qa"
    x, xb, xf = G["x"], G["xb"], G["xf"]
    freq = G["frequencies"]
    control = [analyse(runs[0]["meanMaps"][:, :, k], x) for k in range(len(freq))]
    allrows = []
    summaries = []
    positions = []
    fine = {}
    onsets = []
    for run, ref in [(2, [-0.2, 0.2]), (1, [-0.2, 0])]:
        D = runs[run]
        evidence = []
        for k, f in enumerate(freq):
            A = analyse(D["meanMaps"][:, :, k], x, ref)
            blocks = [analyse(D["macroMaps"][:, :, k, b], xb, ref) for b in range(5)]
            loo = [
                analyse(
                    np.mean(np.delete(D["macroMaps"][:, :, k, :], b, axis=2), axis=2),
                    xb,
                    ref,
                )
                for b in range(5)
            ]
            splits = [
                analyse(np.mean(D["macroMaps"][:, :, k, ids], axis=2), xb, ref)
                for ids in [slice(0, 2), slice(2, 5)]
            ]
            control_conflict = sidelobe_conflict(A, control[k], x[1] - x[0])
            block_ok = [same_pair(A, b) for b in blocks]
            loo_ok = [same_pair(A, b) for b in loo]
            split_ok = [same_pair(A, b) for b in splits]
            # Independent identity-tolerance sensitivity, no dB threshold.
            sensitivity = {
                str(t): all(same_pair(A, b, t) for b in blocks + loo + splits)
                for t in [0.03, 0.04, 0.05]
            }
            temporal = all(block_ok + loo_ok + split_ok)
            robust = A["pair"] is not None and temporal and not control_conflict
            state = 2 if robust else (1 if A["pair"] is not None else 0)
            row = dict(
                run=run,
                separation_m=ref[1] - ref[0],
                frequency_Hz=int(f),
                detected_interior_maxima=sum(not z["boundary"] for z in A["peaks"]),
                mean_pair=A["pair"] is not None,
                saddle_depth_dB=A["saddle_depth_dB"],
                matched_macroblocks=sum(block_ok),
                macroblocks=5,
                leave_one_out_pass=all(loo_ok),
                split_pass=all(split_ok),
                single_source_conflict=control_conflict,
                temporal_support=temporal,
                evidence_stage=state,
                identity_30mm_pass=sensitivity["0.03"],
                identity_50mm_pass=sensitivity["0.05"],
            )
            allrows.append(row)
            evidence.append(row)
            if f in [1000, 3000, 6000]:
                j = int(np.flatnonzero(D["fineF"] == f)[0])
                B = analyse(D["fineMaps"][:, :, j], xf, ref)
                grid = same_pair(A, B) if A["pair"] is not None else B["pair"] is None
                supported = robust and grid
                fine[(run, int(f))] = dict(
                    analysis=B,
                    supported=supported,
                    power=D["fineMaps"][:, :, j],
                    ref=ref,
                )
                s = dict(
                    run=run,
                    separation_m=ref[1] - ref[0],
                    frequency_Hz=int(f),
                    status="ROBUST_DUAL_PEAK" if supported else "UNRESOLVED",
                    saddle_depth_dB=B["saddle_depth_dB"] if supported else None,
                    grid_301_601_consistent=grid,
                    matched_macroblocks=sum(block_ok),
                    reason=(
                        "All deterministic tests pass"
                        if supported
                        else "No pair or insufficient temporal/control/grid evidence"
                    ),
                )
                if supported:
                    pp = D["fineMaps"][:, :, j]
                    s["peak_separation_m"] = float(
                        np.hypot(
                            B["pair"][1]["x_m"] - B["pair"][0]["x_m"],
                            B["pair"][1]["y_m"] - B["pair"][0]["y_m"],
                        )
                    )
                    for q, p in enumerate(B["pair"]):
                        ix = int(np.argmin(abs(xf - p["x_m"])))
                        iy = int(np.argmin(abs(xf - p["y_m"])))
                        s[f"x{q+1}_minima_span_m"] = minima_span(xf, pp[iy, :], ix)
                        s[f"y{q+1}_minima_span_m"] = minima_span(xf, pp[:, ix], iy)
                        err = np.hypot(p["x_m"] - ref[q], p["y_m"])
                        s.update(
                            {
                                f"x{q+1}_m": p["x_m"],
                                f"y{q+1}_m": p["y_m"],
                                f"error{q+1}_m": err,
                            }
                        )
                        positions.append(
                            dict(
                                run=run,
                                separation_m=ref[1] - ref[0],
                                frequency_Hz=int(f),
                                source=q + 1,
                                nominal_x_m=ref[q],
                                nominal_y_m=0,
                                x_m=p["x_m"],
                                y_m=p["y_m"],
                                error_m=err,
                                interpretation="Nominal geometry; physical channel wiring unverified",
                            )
                        )
                summaries.append(s)
        sustained = [
            k
            for k in range(len(evidence) - 2)
            if all(z["evidence_stage"] == 2 for z in evidence[k : k + 3])
        ]
        onset = None if not sustained else int(freq[sustained[0]])
        # Full-grid and finer frequency-neighborhood checks intentionally required before onset claims.
        onsets.append(
            dict(
                run=run,
                separation_m=ref[1] - ref[0],
                f_onset_Hz="indeterminate",
                first_coarse_sustained_candidate_Hz=onset,
                scan_step_Hz=250,
                FFT_bin_Hz=5,
                consecutive_frequencies_required=3,
                later_coarse_losses_Hz=";".join(
                    str(int(z["frequency_Hz"]))
                    for z in evidence
                    if onset and z["frequency_Hz"] > onset and z["evidence_stage"] < 2
                ),
                reason="Coarse candidate only: dense 601-grid neighborhood and time-block independence not established; no physical onset claimed",
            )
        )
    prop = []
    P = runs[3]
    sp = 10 * np.log10(P["spectrumPower"] / (20e-6) ** 2)
    # Reconfirm the existing predeclared selection, not a new acoustic resolution cutoff.
    sf = P["spectrumFrequency"].astype(float)
    peaks = np.flatnonzero((sp[1:-1] > sp[:-2]) & (sp[1:-1] >= sp[2:])) + 1
    prom = np.array(
        [
            sp[i] - np.median(sp[(abs(sf - sf[i]) <= 50) & (abs(sf - sf[i]) > 10)])
            for i in peaks
        ]
    )
    selected = []
    for lo, hi in [(1000, 2000), (2500, 4000), (4500, 6500)]:
        ids = peaks[(sf[peaks] >= lo) & (sf[peaks] <= hi) & (prom >= 6)]
        if not len(ids):
            ids = np.flatnonzero((sf >= lo) & (sf <= hi))
        selected.append(sf[ids[np.argmax(sp[ids])]])
    assert np.array_equal(selected, P["fineF"])
    for k, f in enumerate(P["fineF"]):
        p = P["fineMaps"][:, :, k]
        iy, ix = np.unravel_index(np.argmax(p), p.shape)
        prop.append(
            dict(
                frequency_Hz=int(f),
                SPL_estimate_dB_re20uPa=float(sp[np.argmin(abs(sf - f))]),
                peak_x_m=xf[ix],
                peak_y_m=xf[iy],
                x_minima_span_m=minima_span(xf, p[iy, :], ix),
                y_minima_span_m=minima_span(xf, p[:, ix], iy),
                selection=(
                    "Prominent spectral feature"
                    if k < 2
                    else "Strongest band bin; mechanism unassigned"
                ),
                calibration="Pa metadata; sensitivity calibration chain unverified",
            )
        )
    tables = {
        "T11_resolution_summary": summaries,
        "T11_onset_summary": onsets,
        "T11_position_accuracy": positions,
        "T11_propeller_frequency_summary": prop,
        "T11_frequency_evidence": allrows,
    }
    xy = G["xy"]
    diam = np.max(np.linalg.norm(xy[:, None] - xy[None, :], axis=2))
    tables["T7_acquisition_geometry"] = [
        dict(
            run=r,
            condition=[
                "Single loudspeaker",
                "Two loudspeakers",
                "Two loudspeakers",
                "Two-bladed propeller",
            ][r],
            channels=56,
            sample_count=500000,
            nominal_fs_Hz=50000,
            metadata_fs_Hz=runs[r]["metadata_fs_Hz"],
            duration_s=10,
            source_x_m=["unknown", "-0.20;0.00", "-0.20;0.20", "unknown"][r],
            h_m=1 if r in [1, 2] else "not independently recorded",
            array_Dx_m=np.ptp(xy[:, 0]),
            array_Dy_m=np.ptp(xy[:, 1]),
            array_Dmax_m=diam,
            blades=2 if r == 3 else "",
            rotation_rev_s=100 if r == 3 else "",
            BPF_Hz=200 if r == 3 else "",
            calibration="Unverified chain; stored Pascals",
            channel_order="Sound Pressure_0..55 to coordinate rows 1..56; physical wiring unverified",
        )
        for r in range(4)
    ]
    for name, rows in tables.items():
        csvwrite(out / (name + ".csv"), rows)
    return G, runs, fine, tables


def minima_endpoints(axis, p, index):
    idx = np.flatnonzero((p[1:-1] < p[:-2]) & (p[1:-1] < p[2:])) + 1
    a = idx[idx < index]
    b = idx[idx > index]
    if not len(a) or not len(b):
        return None
    return float(axis[a[-1]]), float(axis[b[0]])


def exact_data(root, geometry, records):
    G, runs, fine, tables = process(root, geometry, records)
    assert G["xy"].shape == (56, 2) and np.unique(G["xy"], axis=0).shape[0] == 56
    assert len(G["xf"]) == 601 and G["xf"][300] == 0
    for row in tables["T11_position_accuracy"]:
        row["e_x_mm"] = 1000 * (row["x_m"] - row["nominal_x_m"])
        row["e_y_mm"] = 1000 * (row["y_m"] - row["nominal_y_m"])
        row["e_mm"] = 1000 * row["error_m"]
        assert np.isclose(
            np.hypot(row["e_x_mm"], row["e_y_mm"]), row["e_mm"], atol=1e-12
        )
    spans = []
    P = runs[3]
    x = G["xf"]
    for k, f in enumerate(P["fineF"]):
        p = P["fineMaps"][:, :, k]
        iy, ix = np.unravel_index(p.argmax(), p.shape)
        for direction, curve, idx in [("x", p[iy, :], ix), ("y", p[:, ix], iy)]:
            endpoints = minima_endpoints(x, curve, idx)
            span = None if endpoints is None else endpoints[1] - endpoints[0]
            approved = tables["T11_propeller_frequency_summary"][k][
                direction + "_minima_span_m"
            ]
            assert (span is None and approved is None) or np.isclose(
                span, approved, atol=1e-14
            )
            spans.append(
                dict(
                    frequency_Hz=int(f),
                    direction=direction,
                    peak_coordinate_m=x[idx],
                    left_minimum_m=None if endpoints is None else endpoints[0],
                    right_minimum_m=None if endpoints is None else endpoints[1],
                    span_m=span,
                    status=(
                        "BOUNDARY_TRUNCATED_OR_UNDEFINED"
                        if span is None
                        else "TWO_FLANKING_MINIMA"
                    ),
                )
            )
    tables["T11_profile_minima_endpoints"] = spans
    for row in tables["T11_frequency_evidence"]:
        f = row["frequency_Hz"]
        key = (row["run"], f)
        row["fine_grid_checked"] = key in fine
        row["fine_grid_pair_supported"] = key in fine and fine[key]["supported"]
    # Subsequent user confirmation supersedes the earlier unknown geometry for 0,3.
    for row in tables["T7_acquisition_geometry"]:
        row["h_m"] = 1.0
        if row["run"] in [0, 3]:
            row["source_x_m"] = "0.00"
        row["source_y_m"] = 0.0
        row["geometry_evidence"] = (
            "Nominal h=1 m; run0,3 on-axis schematic only, not surveyed; dual-source positions from user recollection"
        )
    for name, rows in tables.items():
        csvwrite(root / "figures_tables" / (name + ".csv"), rows)
    return G, runs, fine, tables


# 4. Final figures
BLUE = "#1B619A"
GOLD = "#C89427"
NAVY = "#174875"
GRAY = "#646D74"
GREEN = "#3C7761"
assert "Times New Roman" in {
    f.name for f in font_manager.fontManager.ttflist
}, "Times New Roman required; formal rendering blocked."
font_manager.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "font.size": 10,
        "axes.labelsize": 10.5,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.linewidth": 0.85,
        "lines.linewidth": 1.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "svg.hashsalt": "MECH5961-V5",
        "path.simplify": False,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
        "mathtext.cal": "Times New Roman:italic",
        "mathtext.sf": "Times New Roman",
        "mathtext.tt": "Times New Roman",
        "mathtext.fallback": None,
        "savefig.facecolor": "white",
    }
)
MANIFEST = []
TRACES = {}


def figmm(w, h):
    return plt.figure(figsize=(w / 25.4, h / 25.4))


def tag(fig, ax, text, offset_mm=13):
    b = ax.get_position()
    fig.text(
        (b.x0 + b.x1) / 2,
        b.y0 - offset_mm / (fig.get_figheight() * 25.4),
        text,
        ha="center",
        va="top",
        fontsize=10,
    )


def rel(p):
    assert np.all(p > 0)
    return 10 * np.log10(p / p.max())


def arrow(ax, a, b, color=GRAY):
    ax.annotate(
        "",
        xy=b,
        xytext=a,
        arrowprops=dict(arrowstyle="<->", lw=0.9, color=color, shrinkA=0, shrinkB=0),
    )


def save(fig, name, title, placement, source, method, limits, alignment="equal"):
    fig.canvas.draw()
    boxes = np.array([a.get_position().bounds for a in fig.axes])
    pt = boxes * np.array([fig.get_figwidth() * 72, fig.get_figheight() * 72] * 2)
    if alignment == "equal":
        assert np.ptp(pt[:, 2]) < 1.5 and np.ptp(pt[:, 3]) < 1.5
    elif alignment == "stacked":
        assert np.ptp(pt[:, 0]) < 1.5 and np.ptp(pt[:, 2]) < 1.5
    elif alignment == "mixed":
        assert len(fig.axes) == 2
    record = dict(
        verdict="NOT APPLICABLE" if alignment == "mixed" else "PASS",
        plot_areas_pt=pt.tolist(),
        contract=alignment,
        reason=(
            "Unequal heights serve distinct evidence; F11E xy requires equal physical aspect."
            if alignment != "equal"
            else "Equal comparable areas within 1.5 pt"
        ),
    )
    (ROOT / "qa" / f"{name}.alignment.json").write_text(json.dumps(record, indent=2))
    for ext in ["png", "svg", "pdf"]:
        kw = {"dpi": 600} if ext == "png" else {}
        if ext == "pdf":
            kw["metadata"] = {"CreationDate": None, "ModDate": None}
        if ext == "svg":
            kw["metadata"] = {"Date": None}
        destination = ROOT / "figures_tables" / f"{name}.{ext}"
        temporary = destination.with_name(name + ".staging." + ext)
        fig.savefig(temporary, **kw)
        os.replace(temporary, destination)
    MANIFEST.append(
        dict(
            filename=f"{name}.png;{name}.svg;{name}.pdf",
            english_title=title,
            placement=placement,
            source=source,
            method_parameters_units=method,
            limitations_review=limits,
            width_mm=fig.get_figwidth() * 25.4,
            height_mm=fig.get_figheight() * 25.4,
            dpi=600,
            python_source="postprocessing/run_lab2.py",
            font="Times New Roman",
        )
    )
    plt.close(fig)


def f7a(xy):
    fig = figmm(130, 120)
    ax = fig.add_axes([0.16, 0.16, 0.73, 0.79])
    ax.scatter(xy[:, 0], xy[:, 1], s=17, color=BLUE, zorder=3)
    ax.plot(0, 0, "+", color="#242424", ms=8, mew=1.4)
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    arrow(ax, (xmin, 0.385), (xmax, 0.385))
    ax.text(
        (xmin + xmax) / 2,
        0.4,
        rf"$D_x={xmax-xmin:.3f}$ m",
        ha="center",
        fontsize=10,
        color=GRAY,
    )
    arrow(ax, (0.385, ymin), (0.385, ymax))
    ax.text(
        0.4,
        (ymin + ymax) / 2,
        rf"$D_y={ymax-ymin:.3f}$ m",
        rotation=90,
        rotation_mode="anchor",
        ha="center",
        va="top",
        fontsize=10,
        color=GRAY,
    )
    ax.text(-0.385, -0.38, r"$z=0$", fontsize=10)
    ax.set(
        xlim=(-0.41, 0.45),
        ylim=(-0.41, 0.45),
        xlabel="$x$ (m)",
        ylabel="$y$ (m)",
        xticks=[-0.3, -0.15, 0, 0.15, 0.3],
        yticks=[-0.3, -0.15, 0, 0.15, 0.3],
    )
    ax.set_aspect("equal")
    TRACES["array_xy"] = xy
    save(
        fig,
        "F7A_array_geometry",
        "Measured-coordinate planar microphone geometry",
        "Section 7.1",
        "geometry.mat:xy",
        "All 56 coordinate rows; equal xy aspect; independent coordinate-extrema spans in m",
        "Dx and Dy are not Dmax; physical channel wiring unverified.",
    )


def project(x, y, z):
    return np.asarray(x) + 0.42 * np.asarray(y), 0.38 * np.asarray(
        y
    ) + 0.64 * np.asarray(z)


def f7b(xy):
    fig = figmm(183, 180)
    axes = fig.subplots(2, 2)
    fig.subplots_adjust(
        left=0.045, right=0.96, bottom=0.10, top=0.97, wspace=0.14, hspace=0.35
    )
    for i, ax in enumerate(axes.flat):
        u, v = project(xy[:, 0], xy[:, 1], 0)
        hull = ConvexHull(np.c_[u, v]).vertices
        ax.fill(u[hull], v[hull], fc="#eff3f5", ec="#bdc8d0", lw=0.7)
        ax.scatter(u, v, s=5.5, color=BLUE, zorder=3)
        for target, label, offset in [
            (project(0.48, 0, 0), "x", (0.02, -0.01)),
            (project(0, 0.48, 0), "y", (0.015, 0.015)),
            (project(0, 0, 1.23), "z", (0.02, 0.015)),
        ]:
            ax.annotate(
                "",
                xy=target,
                xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=0.85),
            )
            ax.text(
                target[0] + offset[0],
                target[1] + offset[1],
                f"${label}$",
                fontsize=10,
                color=GRAY,
            )
        ax.text(-0.46, -0.21, r"$z=0$", fontsize=10, color=GRAY)
        if i in [1, 2]:
            refs = [-0.2, 0] if i == 1 else [-0.2, 0.2]
            for j, sx in enumerate(refs):
                px, pv = project(sx, 0, 1)
                ax.scatter(
                    [px],
                    [pv],
                    s=42,
                    marker="o" if j == 0 else "s",
                    color=BLUE if j == 0 else GOLD,
                    zorder=5,
                )
                ax.text(
                    px + (-0.025 if j == 0 else 0.025),
                    pv - 0.10,
                    f"{sx:+.2f} m",
                    ha="right" if j == 0 else "left",
                    fontsize=9,
                    color=BLUE if j == 0 else GOLD,
                )
            arrow(ax, (refs[0], 0.88), (refs[1], 0.88))
            ax.text(
                np.mean(refs),
                0.93,
                rf"$\Delta x={refs[1]-refs[0]:.2f}$ m",
                ha="center",
                fontsize=9,
            )
        else:
            # Requested on-axis schematic; absolute location was not surveyed.
            if i == 0:
                ax.plot(0, 0.64, "D", color=GRAY, ms=6, zorder=5)
            else:
                ax.plot(
                    [-0.065, 0.065],
                    [0.665, 0.615],
                    color=GRAY,
                    lw=3,
                    solid_capstyle="round",
                    zorder=5,
                )
                ax.plot(0, 0.64, "o", color="#252525", ms=3, zorder=6)
        arrow(ax, (-0.50, 0), (-0.50, 0.64))
        ax.text(
            -0.535,
            0.32,
            "$h=1.0$ m",
            rotation=90,
            rotation_mode="anchor",
            ha="center",
            va="center",
            fontsize=9,
        )
        ax.plot([-0.50, -0.24], [0, 0], color="#c6cbd0", lw=0.6)
        ax.plot([-0.50, -0.24], [0.64, 0.64], color="#c6cbd0", lw=0.6)
        ax.set(xlim=(-0.66, 0.61), ylim=(-0.27, 1.07))
        ax.set_aspect("equal")
        ax.axis("off")
        tag(
            fig,
            ax,
            [
                "(a) Single loudspeaker",
                r"(b) $\Delta x=0.20$ m",
                r"(c) $\Delta x=0.40$ m",
                "(d) Two-bladed propeller",
            ][i],
            5,
        )
    save(
        fig,
        "F7B_four_source_geometries",
        "Four source-array configurations",
        "Sections 7.2–7.3",
        "Coor_56Array.txt; nominal dual-source geometry; on-axis single/propeller schematic",
        "Oblique projection u=x+0.42y, v=0.38y+0.64z; h=1 m for all runs; a,d centers shown on z axis schematically",
        "Polygon is coordinate convex hull, not apparatus outline. Source dimensions are symbolic; single/propeller absolute locations are schematic; dual-source nominal positions are recalled, not surveyed.",
    )


def f11a(x, fine):
    fig = figmm(183, 225)
    axes = fig.subplots(3, 2)
    fig.subplots_adjust(
        left=0.105, right=0.975, bottom=0.10, top=0.96, wspace=0.27, hspace=0.65
    )
    for row, f in enumerate([1000, 3000, 6000]):
        for col, run in enumerate([2, 1]):
            ax = axes[row, col]
            d = fine[(run, f)]
            L = rel(d["power"])
            curve = L[300, :]
            ax.plot(x, curve, color=NAVY)
            TRACES[f"dual_run{run}_f{f}"] = np.c_[x, curve]
            for ref in d["ref"]:
                ax.plot([ref, ref], [-25, 1], ls=":", lw=0.65, color="#c4c9ce")
            if d["supported"]:
                a = d["analysis"]
                coords = [p["x_m"] for p in a["pair"]]
                for k, px in enumerate(coords):
                    color = BLUE if k == 0 else GOLD
                    ax.plot(
                        [px, px],
                        [-25, 2.5],
                        ls=(0, (2, 2)),
                        color=color,
                        lw=0.85,
                        alpha=0.7,
                    )
                    ax.plot(px, 2.5, "v", color=color, ms=5, clip_on=False)
                ax.text(
                    -0.45,
                    -30,
                    rf'$\Delta L_s={a["saddle_depth_dB"]:.2f}$ dB',
                    fontsize=9,
                )
                ax.text(
                    -0.45,
                    -36,
                    rf"$\hat x_1={coords[0]:+.3f},\ \hat x_2={coords[1]:+.3f}$ m",
                    fontsize=9,
                )
            ax.set(
                xlim=(-0.5, 0.5),
                ylim=(-40, 4),
                xticks=[-0.4, -0.2, 0, 0.2, 0.4],
                yticks=[-40, -20, 0],
                xlabel="$x$ (m)",
            )
            if col == 0:
                ax.set_ylabel(r"$L_{\mathrm{rel}}$ (dB)")
            tag(
                fig,
                ax,
                f"({chr(97+row*2+col)}) "
                + rf'$\Delta x={d["ref"][1]+.2:.2f}$ m, $f={f/1000:g}$ kHz',
                12,
            )
    save(
        fig,
        "F11A_two_source_profiles",
        "Exact horizontal profiles and two-dimensional peaks",
        "Appendix Q1",
        "run1,2.mat:fineMaps; robust analyzer",
        "601 original samples at y=0; full-2D maximum reference; ΔLs from 8-connected merge saddle",
        "Triangles and colored dashed guides share accepted 2D peak x coordinates; light dotted guides retain nominal x. They are not 1D extrema; y generally differs from zero. No pair/saddle annotations for unresolved conditions.",
    )


def f11b(tables):
    fig = figmm(183, 165)
    ax = fig.add_axes([0.15, 0.56, 0.80, 0.37])
    cat = fig.add_axes([0.15, 0.22, 0.80, 0.17])
    for run, sep, c, m, y in [(2, 0.4, BLUE, "o", 1), (1, 0.2, GOLD, "s", 0)]:
        rows = [r for r in tables["T11_frequency_evidence"] if r["run"] == run]
        for r in rows:
            f = r["frequency_Hz"] / 1000
            s = r["evidence_stage"]
            depth = r["saddle_depth_dB"]
            if depth is not None:
                ax.scatter(
                    f,
                    depth,
                    s=23,
                    marker=m,
                    facecolors=c if s == 2 else "none",
                    edgecolors=c,
                    lw=1,
                    zorder=3,
                )
            if s == 0:
                cat.plot(f, y, "x", color="#a4a8ac", ms=4, mew=1)
            else:
                cat.plot(f, y, m, color=c, mfc=c if s == 2 else "none", ms=5, mew=1)
            if r["fine_grid_pair_supported"]:
                cat.plot(f, y, "o", color="#222222", mfc="none", ms=9, mew=0.75)
        ax.plot([], [], m, color=c, ls="none", label=rf"$\Delta x={sep:.2f}$ m")
    ax.set(
        xlim=(0.85, 6.15),
        ylim=(-0.5, 23),
        xticks=[1, 2, 3, 4, 5, 6],
        yticks=[0, 5, 10, 15, 20],
        ylabel=r"$\Delta L_s$ (dB)",
    )
    ax.legend(loc="upper left", ncol=2, frameon=False)
    tag(fig, ax, "(a) Two-dimensional saddle prominence", 8)
    cat.set(
        xlim=(0.85, 6.15),
        ylim=(-0.55, 1.55),
        xticks=[1, 2, 3, 4, 5, 6],
        yticks=[0, 1],
        yticklabels=[r"$\Delta x=0.20$ m", r"$\Delta x=0.40$ m"],
        xlabel="$f$ (kHz)",
    )
    tag(fig, cat, "(b) Synchronous time-block evidence", 13)
    handles = [
        Line2D([], [], marker="x", ls="none", color="#a4a8ac", label="No pair"),
        Line2D(
            [], [], marker="s", ls="none", color=GRAY, mfc="none", label="Record only"
        ),
        Line2D([], [], marker="s", ls="none", color=GRAY, label="Block-stable"),
        Line2D(
            [],
            [],
            marker="o",
            ls="none",
            color="#222",
            mfc="none",
            ms=8,
            label="601-grid checked",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.53, 0.015),
        ncol=4,
        frameon=False,
        handlelength=1,
        columnspacing=1,
    )
    save(
        fig,
        "F11B_robust_resolvability",
        "Quantitative saddle and discrete stability evidence",
        "Appendix Q1",
        "MicData_0..2.mat: reconstructed mean and synchronized block maps",
        "21 frequencies per spacing, 250 Hz steps; five synchronous 2 s blocks, splits and leave-one-out checks",
        "Undefined saddles absent, not zero. Hollow=record only, filled=block-stable; black ring adds prescribed-frequency 601-grid pair agreement. Most frequencies lack fine-grid check. Onset indeterminate.",
        alignment="stacked",
    )


def f11c(P):
    fig = figmm(183, 105)
    axes = fig.subplots(1, 2)
    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.23, top=0.96, wspace=0.26)
    f = P["spectrumFrequency"].astype(float)
    power = P["spectrumPower"]
    assert np.all(power > 0)
    sp = 10 * np.log10(power / (20e-6) ** 2)
    for j, ax in enumerate(axes):
        use = (f >= (50 if j == 0 else 100)) & (f <= (10000 if j == 0 else 1000))
        scale = 1000 if j == 0 else 1
        ax.plot(f[use] / scale, sp[use], color=NAVY, lw=1.3)
        TRACES["spectrum_" + str(j)] = np.c_[f[use], sp[use]]
        ax.set(ylim=(10, 68), xlabel="$f$ (kHz)" if j == 0 else "$f$ (Hz)")
        if j == 0:
            ax.set(
                xlim=(0.05, 10),
                xticks=[0, 2, 4, 6, 8, 10],
                ylabel="$L_p$ (dB re 20 µPa)",
            )
            for k, hz in enumerate(P["fineF"]):
                c = [GOLD, GRAY, GREEN][k]
                v = sp[np.flatnonzero(f == hz)[0]]
                ax.axvline(
                    hz / 1000, ymax=0.9, lw=0.9, ls=["--", ":", "-."][k], color=c
                )
                ax.plot(hz / 1000, v, ["o", "s", "D"][k], ms=4, color=c)
                ax.text(hz / 1000, 63, f"{hz/1000:g}", ha="center", fontsize=9, color=c)
        else:
            ax.set(xlim=(100, 1000), xticks=[200, 400, 600, 800, 1000])
            for hz, c, label in [
                (200, BLUE, r"$f_{\mathrm{BPF}}$"),
                (400, GOLD, r"$2f_{\mathrm{BPF}}$"),
            ]:
                ax.axvline(
                    hz, ymax=0.82, lw=0.9, ls="--" if hz == 200 else ":", color=c
                )
                ax.text(hz + 15, 61, label, fontsize=10, color=c)
        tag(fig, ax, ["(a) 50 Hz–10 kHz", "(b) 100–1000 Hz"][j], 14)
    save(
        fig,
        "F11C_propeller_spectrum",
        "Microphone 1 spectrum and BPF references",
        "Appendix Q2",
        "MicData_3.mat reconstruction:spectrumFrequency,spectrumPower",
        "All original 5 Hz bins; Hann pressure-bin spectrum; 20 µPa; selected 2000/3000/6185 Hz",
        "Calibration chain unverified. 200/400 Hz guides are theory; anechoic cutoff near 200 Hz. No mechanism assigned to 6185 Hz.",
    )


def f11d(x, P, tables):
    fig = figmm(183, 185)
    axes = fig.subplots(2, 3)
    fig.subplots_adjust(
        left=0.09, right=0.975, bottom=0.12, top=0.97, wspace=0.26, hspace=0.6
    )
    for k, f in enumerate(P["fineF"]):
        p = P["fineMaps"][:, :, k]
        L = rel(p)
        iy, ix = np.unravel_index(p.argmax(), p.shape)
        for row, direction in enumerate(["x", "y"]):
            ax = axes[row, k]
            curve = L[iy, :] if row == 0 else L[:, ix]
            idx = ix if row == 0 else iy
            span = next(
                z
                for z in tables["T11_profile_minima_endpoints"]
                if z["frequency_Hz"] == f and z["direction"] == direction
            )
            ax.plot(x, curve, color=NAVY)
            TRACES[f"prop_{f}_{direction}"] = np.c_[x, curve]
            ax.plot([x[idx], x[idx]], [-25, 0], ls=(0, (2, 2)), lw=0.8, color="#8eaabd")
            ax.plot(x[idx], 0, "o", color=GOLD, ms=4)
            if span["span_m"] is not None:
                a, b = span["left_minimum_m"], span["right_minimum_m"]
                arrow(ax, (a, -29), (b, -29))
                ax.plot([a, a], [-27, -31], color=GRAY, lw=0.8)
                ax.plot([b, b], [-27, -31], color=GRAY, lw=0.8)
            ax.text(
                0.04,
                0.14,
                rf"$\hat {direction}={x[idx]:+.3f}$ m",
                transform=ax.transAxes,
                fontsize=9,
            )
            label = (
                r"$w_{\min}=\mathrm{N/A}$"
                if span["span_m"] is None
                else rf'$w_{{\min}}={span["span_m"]:.3f}$ m'
            )
            ax.text(0.04, 0.04, label, transform=ax.transAxes, fontsize=9, color=GRAY)
            ax.set(
                xlim=(-0.5, 0.5),
                ylim=(-42, 2),
                xticks=[-0.4, 0, 0.4],
                yticks=[-20, -10, 0],
                xlabel=f"${direction}$ (m)",
            )
            if k == 0:
                ax.set_ylabel(r"$L_{\mathrm{rel}}$ (dB)")
            tag(
                fig,
                ax,
                f"({chr(97+row*3+k)}) $f={f/1000:.3f}$ kHz\n{direction}-profile",
                12,
            )
    save(
        fig,
        "F11D_propeller_profiles",
        "Peak-through propeller profiles and lobe spans",
        "Appendix Q2",
        "MicData_3.mat reconstruction:fineMaps",
        "601 samples through dominant xy; strict flanking minima supply actual bracket endpoints",
        "2000 Hz spans undefined; brackets are not centered constructs or half-power widths. Real offsets preserved; width is not rotor diameter.",
    )


def f11e(rows):
    fig = figmm(183, 125)
    left = fig.add_axes([0.09, 0.35, 0.43, 0.44])
    right = fig.add_axes([0.64, 0.24, 0.33, 0.62])
    for r in rows:
        c = BLUE if r["source"] == 1 else GOLD
        m = "o" if r["source"] == 1 else "s"
        a, b = [r["nominal_x_m"], r["x_m"]], [r["nominal_y_m"], r["y_m"]]
        left.plot(a, b, color=c, lw=1)
        left.plot(a[0], b[0], m, mfc="white", mec=c, ms=6, mew=1.1)
        left.plot(a[1], b[1], m, color=c, ms=4)
        group = {(2, 3000): 0, (2, 6000): 1, (1, 6000): 2}[
            (r["run"], r["frequency_Hz"])
        ]
        for key, mk, size, face in [
            ("e_x_mm", "o", 4, c),
            ("e_y_mm", "s", 7, "none"),
            ("e_mm", "^", 5, c),
        ]:
            right.plot(group, r[key], mk, ms=size, color=c, mfc=face, mew=1)
    left.set(
        xlim=(-0.27, 0.28),
        ylim=(-0.045, 0.10),
        xticks=[-0.2, 0, 0.2],
        yticks=[0, 0.04, 0.08],
        xlabel="$x$ (m)",
        ylabel="$y$ (m)",
    )
    left.set_aspect("equal")
    left.legend(
        handles=[
            Line2D(
                [], [], marker="o", color=GRAY, ls="none", mfc="white", label="Nominal"
            ),
            Line2D([], [], marker="o", color=GRAY, ls="none", label="Identified"),
        ],
        loc="upper left",
        bbox_to_anchor=(0, 1.7),
        frameon=False,
        ncol=2,
    )
    right.axhline(0, color="#acb0b3", lw=0.8)
    right.set(
        xlim=(-0.4, 2.4),
        ylim=(-8, 48),
        xticks=[0, 1, 2],
        xticklabels=["0.40 m\n3 kHz", "0.40 m\n6 kHz", "0.20 m\n6 kHz"],
        ylabel="Nominal residual (mm)",
    )
    right.legend(
        handles=[
            Line2D([], [], marker=m, color=GRAY, ls="none", mfc=face, label=lab)
            for m, face, lab in [
                ("o", GRAY, "$e_x$"),
                ("s", "white", "$e_y$"),
                ("^", GRAY, "$e$"),
            ]
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=3,
        frameon=False,
        columnspacing=0.7,
        handlelength=1,
    )
    tag(fig, left, "(a) Nominal and identified positions", 14)
    tag(fig, right, "(b) Signed components and magnitude", 18)
    fig.legend(
        handles=[
            Line2D([], [], marker="o", color=BLUE, ls="none", label="Source 1"),
            Line2D([], [], marker="s", color=GOLD, ls="none", label="Source 2"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=2,
        frameon=False,
    )
    save(
        fig,
        "F11E_nominal_error",
        "Nominal xy positions and signed residuals",
        "Appendix Q3 — review candidate",
        "Recomputed T11_position_accuracy.csv",
        "Equal-aspect xy; ex=xhat-xnom, ey=yhat-ynom, e=hypot(ex,ey); six residuals in mm",
        "Only three supported conditions; nominal geometry and unverified wiring prevent certified accuracy. Coincident markers retained, not jittered.",
        alignment="mixed",
    )


def render(geometry, records):
    G, runs, fine, tables = exact_data(ROOT, geometry, records)
    f7a(G["xy"])
    f7b(G["xy"])
    f11a(G["xf"], fine)
    f11b(tables)
    f11c(runs[3])
    f11d(G["xf"], runs[3], tables)
    f11e(tables["T11_position_accuracy"])
    for name in tables:
        MANIFEST.append(
            dict(
                filename=name + ".csv",
                english_title=name.replace("_", " "),
                placement="Section 7" if name.startswith("T7") else "Appendix Q1–Q3",
                source="Current-run original MAT reconstruction",
                method_parameters_units="Recomputed from original pressure; units in columns",
                limitations_review="Undefined metrics empty; one record per condition; calibration and wiring unverified",
            )
        )
    csvwrite(ROOT / "figures_tables" / "FIGURE_TABLE_MANIFEST.csv", MANIFEST)
    np.savez_compressed(ROOT / "qa" / "V5_plotted_samples.npz", **TRACES)
    make_overview(ROOT, MANIFEST[:7])
    print("Seven final figures regenerated.")
    return tables


def make_overview(root, manifest):
    fig, axs = plt.subplots(4, 2, figsize=(12, 16))
    fig.subplots_adjust(wspace=0.02, hspace=0.08)
    for ax, item in zip(axs.flat, manifest):
        ax.imshow(plt.imread(root / "figures_tables" / item["filename"].split(";")[0]))
        ax.axis("off")
    axs.flat[-1].axis("off")
    fig.savefig(root / "qa" / "figure_overview.png", dpi=150)
    plt.close(fig)


# 5. Verification and transactional publication
def verify_outputs(root, tables):
    expected = {
        (2, 1000): "UNRESOLVED",
        (2, 3000): "ROBUST_DUAL_PEAK",
        (2, 6000): "ROBUST_DUAL_PEAK",
        (1, 1000): "UNRESOLVED",
        (1, 3000): "UNRESOLVED",
        (1, 6000): "ROBUST_DUAL_PEAK",
    }
    assert {
        (r["run"], r["frequency_Hz"]): r["status"]
        for r in tables["T11_resolution_summary"]
    } == expected
    assert all(r["f_onset_Hz"] == "indeterminate" for r in tables["T11_onset_summary"])
    for extension in ["png", "svg", "pdf"]:
        assert len(list((root / "figures_tables").glob("*." + extension))) == 7
    for pdf in (root / "figures_tables").glob("*.pdf"):
        with pymupdf.open(pdf) as document:
            for page in document:
                for font in page.get_fonts():
                    assert "TimesNewRoman" in font[3].replace(" ", ""), (pdf, font)
                    assert document.extract_font(font[0])[3], (pdf, "font not embedded")


def input_hashes(lab_root):
    inputs = list((lab_root / "Beamforming").rglob("*")) + list(
        (lab_root / "Group2-20261116T100220Z-1-001").rglob("*")
    )
    result = {}
    for path in inputs:
        if path.is_file():
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                    digest.update(chunk)
            result[path.relative_to(lab_root).as_posix()] = digest.hexdigest()
    return result


def main():
    global ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, help="Optional alternative final output directory."
    )
    arguments = parser.parse_args()
    lab_root = Path(__file__).resolve().parent.parent
    destination = (
        arguments.output.resolve() if arguments.output else lab_root / "figures_tables"
    )
    required = [lab_root / "Beamforming/Coor_56Array.txt"] + [
        lab_root / f"Group2-20261116T100220Z-1-001/MicData_{r}.mat" for r in range(4)
    ]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(
                f"Required external experimental input missing: {path}"
            )
    original_hashes = input_hashes(lab_root)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".lab2-runtime-", dir=destination.parent
    ) as directory:
        ROOT = Path(directory)
        (ROOT / "figures_tables").mkdir()
        (ROOT / "qa").mkdir()
        geometry, records = reconstruct(lab_root)
        tables = render(geometry, records)
        verify_outputs(ROOT, tables)
        assert original_hashes == input_hashes(
            lab_root
        ), "Input checksum changed during execution"
        # Only these script-owned files may be replaced; unknown user files are protected.
        generated = {p.name for p in (ROOT / "figures_tables").iterdir()}
        if destination.exists():
            unknown = {p.name for p in destination.iterdir()} - generated
            if unknown:
                raise RuntimeError(
                    f"Output directory contains unowned files: {sorted(unknown)}. Select --output."
                )
        backup = Path(directory) / "previous-output"
        if destination.exists():
            destination.rename(backup)
        try:
            (ROOT / "figures_tables").rename(destination)
        except Exception:
            if backup.exists():
                backup.rename(destination)
            raise
        print(
            f"PASS: original inputs unchanged; six resolution states verified; seven figures and formal CSVs: {destination}"
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAILED: {type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(1)
