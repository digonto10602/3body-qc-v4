#!/usr/bin/env python3
"""
L=20 KKpi spectrum comparison plotter with recent non-interacting levels.

Panels: 000_A1m, 100_A2, 110_A2, 111_A2, 200_A2
  x=18 : old QC fit spectrum loaded from the old QC_states_jackknifed path
  x=20 : lattice spectrum V20_<irrep>, labelled L=20
  x=22 : fitted model true-zero spectrum from model_true_zeros_all_windows.csv,
         plotted dark orange when --model-true-zeros-file is provided.

If no model/fitted-spectrum file is provided, the script still plots old QC,
lattice, thresholds/cutoff, and recent non-interacting curves without model points.

Non-interacting levels are generated for Lbyas = 16,17,...,24 using the recent
nonint3body_general_doublecover.py software, then plotted as black curves.
The x-axis frame is [16,24], but no tick/label is drawn at x=16 or x=24.
"""

import argparse
import csv
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from sys import platform

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# Thresholds copied from the user's original plotting script.
KKPI = 0.26302
KKPIPI = 0.33208

# Paths copied from the user's original plotting script.
threebody_path_ubuntu = '/home/digonto/Codes/Practical_Lattice_v2/3body_quantization_v2/'
threebody_path_macos = '../../../'

if platform == "linux" or platform == "linux2":
    threebody_path = threebody_path_ubuntu
elif platform == "darwin":
    threebody_path = threebody_path_macos
else:
    threebody_path = threebody_path_ubuntu

IRREPS = ['000_A1m', '100_A2', '110_A2', '111_A2', '200_A2']
PMOMS  = ['000',     '100',    '110',   '111',   '200']
P_TO_OUTPUT = {'000': 'P000', '100': 'P001', '110': 'P110', '111': 'P111', '200': 'P200'}
IRREP_TO_SOFTWARE = {'000_A1m': 'A1m', '100_A2': 'A2', '110_A2': 'A2', '111_A2': 'A2', '200_A2': 'A2'}

# Fixed y-axis ticks/labels requested by user; these are placed on the left axis.
YTICKS = [0.25, KKPI, 0.28, 0.31, KKPIPI, 0.34]
YTICK_LABELS = ['0.25', r'$KK\pi$', '0.28', '0.31', r'$KK\pi\pi$', '0.34']


def safe_usetex(enable=True):
    plt.rcParams.update({'font.size': 22, 'legend.fontsize': 12})
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
    plt.rc('text', usetex=bool(enable))


def resolve_file(path, description, allow_missing=False):
    if path and os.path.exists(path):
        return path
    if allow_missing:
        print(f"[warn] missing optional {description}: {path}")
        return None
    raise FileNotFoundError(f"Could not find {description}: {path}")


def resolve_recent_fit_file(user_path=None):
    if user_path and os.path.exists(user_path):
        return user_path
    patterns = [
        "output_v32w_QC_fitter_norm_refine_v2_fixedcoarse/*bestfit_QC_spectrum.dat",
        "output_v32w*/**/*bestfit_QC_spectrum.dat",
        "**/*bestfit_QC_spectrum.dat",
    ]
    candidates = []
    for pat in patterns:
        candidates.extend(glob.glob(pat, recursive=True))
    candidates = sorted(set(candidates))
    if candidates:
        v32w = [c for c in candidates if 'v32w' in c]
        return v32w[0] if v32w else candidates[0]
    if user_path:
        raise FileNotFoundError(f"Could not find recent fit spectrum file: {user_path}")
    raise FileNotFoundError("Could not auto-find '*bestfit_QC_spectrum.dat'. Pass --recent-fit-file.")


def find_nonint_script(user_path=None):
    if user_path and os.path.exists(user_path):
        return user_path
    candidates = [
        "scripts/nonint3body_general_doublecover.py",
        "source/nonint3body_general_doublecover.py",
        "../nonint3body_general_doublecover_package/source/nonint3body_general_doublecover.py",
        "../../nonint3body_general_doublecover_package/source/nonint3body_general_doublecover.py",
        os.path.expanduser("~/Codes/Practical_Lattice_v2/3body_quantization_v3/studies/nonint3body_general_doublecover_package/source/nonint3body_general_doublecover.py"),
        os.path.expanduser("~/Codes/Practical_Lattice_v2/3body_quantization_v3/nonint3body_general_doublecover_package/source/nonint3body_general_doublecover.py"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    raise FileNotFoundError(
        "Could not find nonint3body_general_doublecover.py. Pass --nonint-script /path/to/source/nonint3body_general_doublecover.py."
    )


def read_old_qc_file(cutoff_val, pmom):
    qc_drive = (
        threebody_path
        + "/test_files/QC_states_twoparams_fitting_L_20_24_v2/cutoff_"
        + str(cutoff_val)
        + "/"
    )
    qc_file = (
        qc_drive
        + "QC_states_jackknifed_L_20_nP_"
        + pmom
        + "_energycutoff_"
        + str(cutoff_val)
        + "_two_params.dat"
    )
    qc_file = resolve_file(qc_file, f"old QC file for nP={pmom}")
    data = np.genfromtxt(qc_file, unpack=True)
    if data.ndim == 1:
        data = data.reshape((len(data), 1))
    qc_spec = np.asarray(data[1], dtype=float)
    qc_err = np.asarray(data[2], dtype=float) if data.shape[0] > 2 else np.zeros_like(qc_spec)
    return qc_spec, qc_err, qc_file


def read_lattice_v20(irrep):
    spec_drive = threebody_path + "/lattice_data/KKpi_interacting_spectrum/twoptvar_analysis/masses/"
    spectrum_file = spec_drive + "V20_" + str(irrep)
    spectrum_file = resolve_file(spectrum_file, f"V20 lattice file for {irrep}")
    data = np.genfromtxt(spectrum_file, unpack=True)
    if data.ndim == 1:
        data = data.reshape((len(data), 1))
    L = np.asarray(data[0], dtype=float)
    E = np.asarray(data[1], dtype=float)
    Eerr = np.asarray(data[3], dtype=float) if data.shape[0] > 3 else np.asarray(data[2], dtype=float)
    return L, E, Eerr, spectrum_file


def read_recent_bestfit_spectrum(path, true_zero_only=True):
    """
    Parse v32w bestfit_QC_spectrum.dat.
    Expected columns:
      label index Ecm kind reason BL BR FL FR yBL_scaled yBR_scaled yFL_scaled yFR_scaled
    """
    out = {ir: [] for ir in IRREPS}
    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            label = parts[0]
            if label not in out:
                continue
            try:
                ecm = float(parts[2])
            except ValueError:
                continue
            kind = parts[3]
            if true_zero_only and kind != 'true_zero':
                continue
            out[label].append(ecm)
    for key in out:
        out[key] = np.array(sorted(out[key]), dtype=float)
    return out


def _parse_bool(value):
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n"}:
        return False
    return None


def _parse_float(value, default=None):
    try:
        if value is None or str(value).strip() == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def read_model_true_zeros_csv(path, fit_cutoff=None, lbyas=20, emax=0.36):
    """
    Parse the new fitter output model_true_zeros_all_windows.csv.

    Expected useful columns include:
      fit_cutoff, true_zero_search_Ecm_max, Lbyas, irrep, model_Ecm,
      user_label, status, used_in_fit_chi2, inside_fit_cutoff, bracket_id, ...

    If the file contains multiple fit_cutoff values and --fit-cutoff is not set,
    select the largest fit_cutoff present. This matches the usual "best/highest
    completed cutoff" plotting workflow while still printing the selected value.
    """
    path = resolve_file(path, "model_true_zeros_all_windows.csv")
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"No CSV header found in {path}")
        fieldnames = [name.strip() for name in reader.fieldnames]
        for raw in reader:
            # Normalize keys while preserving values.
            row = {str(k).strip(): v for k, v in raw.items() if k is not None}
            rows.append(row)

    if not rows:
        return {ir: np.array([], dtype=float) for ir in IRREPS}, None, path

    # Determine available fit cutoffs.
    available_cutoffs = sorted({
        c for c in (_parse_float(r.get("fit_cutoff")) for r in rows)
        if c is not None
    })
    selected_cutoff = fit_cutoff
    if selected_cutoff is None and available_cutoffs:
        selected_cutoff = available_cutoffs[-1]
        print(f"[model] --fit-cutoff not provided; using highest fit_cutoff={selected_cutoff:g} from {path}")
    elif selected_cutoff is not None:
        print(f"[model] using requested fit_cutoff={selected_cutoff:g} from {path}")

    out = {ir: [] for ir in IRREPS}
    kept = 0
    skipped = 0
    for r in rows:
        irrep = (r.get("irrep") or r.get("label") or "").strip()
        if irrep not in out:
            skipped += 1
            continue

        row_L = _parse_float(r.get("Lbyas"), default=float(lbyas))
        if row_L is not None and int(round(row_L)) != int(lbyas):
            skipped += 1
            continue

        row_cut = _parse_float(r.get("fit_cutoff"))
        if selected_cutoff is not None and row_cut is not None and abs(row_cut - selected_cutoff) > 5e-7:
            skipped += 1
            continue

        user_label = _parse_bool(r.get("user_label"))
        if user_label is False:
            skipped += 1
            continue

        status = str(r.get("status") or "").strip().upper()
        if "FAILED" in status or "FALSE" == status:
            skipped += 1
            continue

        ecm = _parse_float(r.get("model_Ecm"))
        if ecm is None:
            skipped += 1
            continue
        if ecm > emax + 1e-12:
            skipped += 1
            continue

        out[irrep].append(ecm)
        kept += 1

    for key in out:
        out[key] = np.array(sorted(out[key]), dtype=float)

    print(f"[model] model_true_zeros_file = {path}")
    print(f"[model] selected_fit_cutoff = {selected_cutoff}")
    print(f"[model] true_zero_search_Ecm_max = {emax}")
    print(f"[model] kept_rows = {kept}; skipped_rows = {skipped}")
    return out, selected_cutoff, path


def write_nonint_config(path, outdir, tag, Lbyas, Emin, Emax, xi, nmax, mK, mpi):
    with open(path, 'w') as f:
        f.write(f"Ecm_min = {Emin}\n")
        f.write(f"Ecm_max = {Emax}\n")
        f.write(f"L = {Lbyas}\n")
        f.write(f"xi = {xi}\n")
        f.write(f"nmax = {nmax}\n")
        f.write("energy_group_tol = 1e-10\n")
        f.write("irrep_tol = 1e-8\n")
        f.write("double_cover = y\n")
        f.write("print_empty_irreps = n\n")
        f.write("P_list = 000 001 110 111 200\n")
        f.write(f"particle1_name = K\nparticle1_mass = {mK}\nparticle1_J = 0\nparticle1_P = -\nparticle1_C = none\n")
        f.write(f"particle2_name = K\nparticle2_mass = {mK}\nparticle2_J = 0\nparticle2_P = -\nparticle2_C = none\n")
        f.write(f"particle3_name = pi\nparticle3_mass = {mpi}\nparticle3_J = 0\nparticle3_P = -\nparticle3_C = -\n")
        f.write(f"output_dir = {outdir}\n")
        f.write(f"output_tag = {tag}\n")


def parse_nonint_output_file(path):
    levels = []
    if not os.path.exists(path):
        return levels
    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                levels.append(float(parts[1]))
            except ValueError:
                continue
    return sorted(levels)


def generate_recent_nonint_levels(args):
    """
    Use the recent nonint3body_general_doublecover.py software for Lbyas=Lmin..Lmax.
    Returns dict: irrep -> list of (L, Ecm) points. Each level index is connected
    across L by sorting levels at each L and using the rank index.
    """
    nonint_script = find_nonint_script(args.nonint_script)
    cache_dir = args.nonint_cache_dir or os.path.join(args.output_nonint_dir, "recent_nonint_L16to24_cache")
    os.makedirs(cache_dir, exist_ok=True)
    print(f"[nonint] using recent software: {nonint_script}")
    print(f"[nonint] cache/output dir: {cache_dir}")

    levels_by_irrep_by_L = {ir: {} for ir in IRREPS}
    for Lbyas in range(args.Lmin, args.Lmax + 1):
        tag = f"recent_nonint_L{Lbyas}"
        config_path = os.path.join(cache_dir, f"{tag}.in")
        marker_path = os.path.join(cache_dir, f"{tag}.done")
        write_nonint_config(
            config_path, cache_dir, tag, Lbyas,
            args.ymin, args.ymax, args.xi, args.nmax, args.mK, args.mpi
        )
        if not os.path.exists(marker_path) or args.force_recompute_nonint:
            print(f"[nonint] running Lbyas={Lbyas}")
            cmd = [sys.executable, nonint_script, config_path]
            subprocess.run(cmd, check=True)
            with open(marker_path, 'w') as f:
                f.write("done\n")
        else:
            print(f"[nonint] reuse Lbyas={Lbyas}")

        for irrep, pmom in zip(IRREPS, PMOMS):
            ptag = P_TO_OUTPUT[pmom]
            soft_ir = IRREP_TO_SOFTWARE[irrep]
            path = os.path.join(cache_dir, f"{tag}_{ptag}_{soft_ir}_nonint3body.dat")
            levels = parse_nonint_output_file(path)
            levels_by_irrep_by_L[irrep][float(Lbyas)] = levels
            print(f"[nonint] L={Lbyas} {irrep}: levels={len(levels)} file={os.path.basename(path)}")

    # Convert sorted levels by L into continuous rank curves.
    curves = {ir: [] for ir in IRREPS}
    Lvals = [float(L) for L in range(args.Lmin, args.Lmax + 1)]
    for irrep in IRREPS:
        max_n = max((len(levels_by_irrep_by_L[irrep].get(L, [])) for L in Lvals), default=0)
        for idx in range(max_n):
            xs, ys = [], []
            for L in Lvals:
                levels = levels_by_irrep_by_L[irrep].get(L, [])
                if idx < len(levels):
                    xs.append(L)
                    ys.append(levels[idx])
            if len(xs) >= 2:
                curves[irrep].append((np.array(xs), np.array(ys)))
    return curves


def plot_spectrum(args):
    safe_usetex(not args.no_tex)
    os.makedirs(args.output_nonint_dir, exist_ok=True)

    selected_fit_cutoff = None
    model_source_file = None
    if args.model_true_zeros_file:
        model_spectrum, selected_fit_cutoff, model_source_file = read_model_true_zeros_csv(
            args.model_true_zeros_file,
            fit_cutoff=args.fit_cutoff,
            lbyas=args.Lmodel,
            emax=args.model_emax,
        )
    elif args.recent_fit_file:
        # Legacy support for old *_bestfit_QC_spectrum.dat files.
        model_source_file = resolve_recent_fit_file(args.recent_fit_file)
        model_spectrum = read_recent_bestfit_spectrum(
            model_source_file, true_zero_only=not args.include_nonzero_recent
        )
    else:
        print("[model] no --model-true-zeros-file or --recent-fit-file provided; plotting without fitted model points")
        model_spectrum = {ir: np.array([], dtype=float) for ir in IRREPS}

    nonint_curves = generate_recent_nonint_levels(args)

    fig, ax = plt.subplots(1, 5, figsize=(24, 10), sharey=False)

    print(f"[input] threebody_path = {threebody_path}")
    print(f"[input] model source file = {model_source_file}")
    print(f"[input] selected fit cutoff = {selected_fit_cutoff}")
    print(f"[config] cutoff guide = {args.cutoff}")

    for i, (irrep, pmom) in enumerate(zip(IRREPS, PMOMS)):
        old_E, old_err, old_file = read_old_qc_file(args.cutoff, pmom)
        L, lat_E, lat_err, lat_file = read_lattice_v20(irrep)
        model_E = model_spectrum.get(irrep, np.array([], dtype=float))

        print(f"[panel {irrep}] old_QC={len(old_E)} lattice={len(lat_E)} model_true_zero={len(model_E)} nonint_curves={len(nonint_curves.get(irrep, []))}")
        print(f"  old_file     = {old_file}")
        print(f"  lattice_file = {lat_file}")

        a = ax[i]
        a.set_title(irrep, fontsize=30)
        a.set_xlim(float(args.Lmin), float(args.Lmax))
        # Hard-fix the y frame and disable autoscale on y so later artists cannot change it.
        a.set_ylim(args.ymin, args.ymax)
        a.set_autoscaley_on(False)
        # User request: frame begins at 16 and ends at 24, but no tick/label at 16 or 24.
        a.set_xticks([18, 20, 22])
        a.set_xticklabels(['old QC', r'$L=20$', 'model fit'], rotation=45, ha='right')

        # Recent non-interacting levels calculated for L=16..24.
        for xs, ys in nonint_curves.get(irrep, []):
            a.plot(xs, ys, linestyle='-', color='black', linewidth=1.5, zorder=2)

        # Threshold and cutoff reference lines.
        a.axhline(y=KKPI, color='grey', linestyle='--', linewidth=1, zorder=1)
        a.axhline(y=KKPIPI, color='grey', linestyle='--', linewidth=1, zorder=1)
        a.axhline(y=args.cutoff, linestyle='--', color='darkred', linewidth=2.2, zorder=3)

        # Left-axis tick marks and labels exactly as requested.
        # Use FixedLocator/FixedFormatter so matplotlib cannot replace them with auto ticks.
        a.yaxis.set_major_locator(mticker.FixedLocator(YTICKS))
        if i == 0:
            a.yaxis.set_major_formatter(mticker.FixedFormatter(YTICK_LABELS))
            a.set_ylabel(r'$a_t E_{\mathrm{cm}}$', fontsize=30)
            a.tick_params(axis='y', which='major', labelleft=True, labelsize=26)
        else:
            a.yaxis.set_major_formatter(mticker.NullFormatter())
            a.tick_params(axis='y', which='major', labelleft=False)
        a.yaxis.set_minor_locator(mticker.NullLocator())

        # Old QC fit spectrum at x=18.
        if len(old_E) > 0:
            x_old = np.full(len(old_E), 18.0)
            if args.show_errors and np.any(old_err > 0):
                a.errorbar(
                    x_old, old_E, yerr=old_err, fmt='o', markerfacecolor='white',
                    markersize=13, color='grey', markeredgewidth=2, capsize=8,
                    linestyle='none', zorder=6,
                )
            else:
                a.plot(
                    x_old, old_E, marker='o', markerfacecolor='white', markersize=13,
                    color='grey', linestyle='none', markeredgewidth=2, zorder=6,
                )

        # Lattice values at x=20.
        x_lat = np.full(len(lat_E), 20.0)
        a.errorbar(
            x_lat, lat_E, yerr=lat_err, fmt='o', markerfacecolor='white',
            markersize=15, color='black', markeredgewidth=2, capsize=10,
            linestyle='none', zorder=7,
        )

        # Fitted model true-zero spectrum at x=22, all levels dark orange.
        if len(model_E) > 0:
            x_model = np.full(len(model_E), 22.0)
            a.plot(
                x_model, model_E, marker='o', markerfacecolor='darkorange',
                markersize=13, color='darkorange', linestyle='none',
                markeredgewidth=1.8, zorder=8,
            )

        a.tick_params(axis='x', which='major', labelsize=18)

    # Final pass: hard-enforce identical y-limits and fixed custom left ticks after all plotting.
    for i, a in enumerate(ax):
        a.set_ylim(args.ymin, args.ymax)
        a.set_autoscaley_on(False)
        a.yaxis.set_major_locator(mticker.FixedLocator(YTICKS))
        if i == 0:
            a.yaxis.set_major_formatter(mticker.FixedFormatter(YTICK_LABELS))
            a.tick_params(axis='y', which='major', labelleft=True, labelsize=26)
        else:
            a.yaxis.set_major_formatter(mticker.NullFormatter())
            a.tick_params(axis='y', which='major', labelleft=False)
        a.yaxis.set_minor_locator(mticker.NullLocator())

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    plt.savefig(args.output, dpi=args.dpi)
    if args.output_pdf:
        plt.savefig(args.output_pdf)
    if args.show:
        plt.show()
    plt.close()
    print(f"[done] saved {args.output}")
    if args.output_pdf:
        print(f"[done] saved {args.output_pdf}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cutoff', type=float, default=0.335,
                    help='Fit cutoff shown as horizontal dark-red dashed line and used in old QC path.')
    ap.add_argument('--model-true-zeros-file', default=None,
                    help='Path to model_true_zeros_all_windows.csv or combined all-cutoffs CSV from the new fitter. If omitted, no model-fit points are plotted unless --recent-fit-file is supplied.')
    ap.add_argument('--fit-cutoff', type=float, default=None,
                    help='When --model-true-zeros-file contains multiple fit_cutoff values, select this one. Default: highest fit_cutoff in the file.')
    ap.add_argument('--model-emax', type=float, default=0.36,
                    help='Maximum model_Ecm to plot from model_true_zeros_all_windows.csv. Default: 0.36.')
    ap.add_argument('--Lmodel', type=int, default=20,
                    help='Lbyas value to select from model_true_zeros_all_windows.csv. Default: 20.')
    ap.add_argument('--recent-fit-file', default=None,
                    help='Optional legacy path to v32w *_bestfit_QC_spectrum.dat. Used only if --model-true-zeros-file is omitted.')
    ap.add_argument('--nonint-script', default=None,
                    help='Path to recent source/nonint3body_general_doublecover.py. Required if auto-detect fails.')
    ap.add_argument('--nonint-cache-dir', default=None,
                    help='Directory for generated recent non-interacting files.')
    ap.add_argument('--output-nonint-dir', default='output_recent_nonint_plotter')
    ap.add_argument('--force-recompute-nonint', action='store_true')
    ap.add_argument('--Lmin', type=int, default=16)
    ap.add_argument('--Lmax', type=int, default=24)
    ap.add_argument('--mK', type=float, default=0.09698)
    ap.add_argument('--mpi', type=float, default=0.06906)
    ap.add_argument('--xi', type=float, default=3.444)
    ap.add_argument('--nmax', type=int, default=5)
    ap.add_argument('--output', default='KKpi_L16to24_oldQC_L20_modelTrueZeros_recentNonint_leftaxis.png')
    ap.add_argument('--output-pdf', default='KKpi_L16to24_oldQC_L20_modelTrueZeros_recentNonint_leftaxis.pdf')
    ap.add_argument('--ymin', type=float, default=0.25)
    ap.add_argument('--ymax', type=float, default=0.35,
                    help='Upper y-limit. Default fixed to 0.35 as requested.')
    ap.add_argument('--dpi', type=int, default=300)
    ap.add_argument('--show', action='store_true')
    ap.add_argument('--no-tex', action='store_true')
    ap.add_argument('--show-errors', action='store_true', help='Also show errors on old QC points if present.')
    ap.add_argument('--include-nonzero-recent', action='store_true', help='Plot all rows in recent bestfit file, not only kind=true_zero.')
    args = ap.parse_args()
    plot_spectrum(args)


if __name__ == '__main__':
    main()
