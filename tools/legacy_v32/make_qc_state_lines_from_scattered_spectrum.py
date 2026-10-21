#!/usr/bin/env python3
"""
make_qc_state_lines_from_scattered_spectrum.py

Read scattered QC spectrum files from v32x/v32y output directories or zip files,
track levels in each irrep as approximately continuous branches in Lbyas, and
write dense interpolated state-line files that can be plotted like non-interacting
spectrum lines.

Supported input files:
  - *_LbyasXX_QC_spectrum.dat with columns:
      Lbyas spatial_L label state_index Ecm kind reason ...
  - *_spectrum_vs_L_all_levels.dat with columns:
      Lbyas spatial_L label state_index Ecm

The script ignores *_all_candidates_*.dat by default.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
import tempfile
import zipfile
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

try:
    from scipy.interpolate import CubicSpline, PchipInterpolator, interp1d
    from scipy.optimize import linear_sum_assignment
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False
    CubicSpline = PchipInterpolator = interp1d = None
    linear_sum_assignment = None


@dataclass
class Point:
    Lbyas: float
    spatial_L: float
    label: str
    Ecm: float
    state_index: int = -1
    source: str = ""
    kind: str = "true_zero"


@dataclass
class Track:
    label: str
    track_id: int
    points: List[Point] = field(default_factory=list)

    def last(self) -> Point:
        return self.points[-1]

    def add(self, p: Point) -> None:
        self.points.append(p)

    def sorted_points(self) -> List[Point]:
        return sorted(self.points, key=lambda p: (p.Lbyas, p.Ecm))


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Track scattered QC spectrum points versus Lbyas and write continuous interpolated state-line files."
    )
    ap.add_argument(
        "inputs",
        nargs="+",
        help="Input output directories, .dat files, or .zip files from v32x/v32y.",
    )
    ap.add_argument("--outdir", default="qc_state_lines_continuous", help="Output directory.")
    ap.add_argument("--lmin", type=float, default=16.0, help="Minimum Lbyas to use/output.")
    ap.add_argument("--lmax", type=float, default=24.0, help="Maximum Lbyas to use/output.")
    ap.add_argument("--xi", type=float, default=3.444, help="Anisotropy xi, spatial_L = Lbyas * xi.")
    ap.add_argument(
        "--irreps",
        nargs="*",
        default=["000_A1m", "100_A2", "110_A2", "111_A2", "200_A2"],
        help="Irreps to include.",
    )
    ap.add_argument(
        "--interp",
        choices=["linear", "cubic", "pchip"],
        default="pchip",
        help="Interpolation method for dense lines. pchip is cubic-like but shape-preserving.",
    )
    ap.add_argument("--dense-n", type=int, default=401, help="Dense output points per state branch.")
    ap.add_argument(
        "--max-jump",
        type=float,
        default=0.035,
        help="Maximum |Delta Ecm| allowed when connecting a state between neighboring integer L values.",
    )
    ap.add_argument(
        "--max-gap",
        type=float,
        default=1.5,
        help="Maximum Lbyas gap allowed for branch matching. Default connects adjacent integer L only.",
    )
    ap.add_argument(
        "--min-points",
        type=int,
        default=2,
        help="Minimum observed points needed to write an interpolated state line.",
    )
    ap.add_argument(
        "--require-full-span",
        action="store_true",
        help="Only write branches that have support at both lmin and lmax.",
    )
    ap.add_argument(
        "--keep-above-ymax",
        action="store_true",
        help="Do not filter by ymax. By default no ymax is applied unless --ymax is given.",
    )
    ap.add_argument("--ymin", type=float, default=None, help="Optional minimum Ecm filter.")
    ap.add_argument("--ymax", type=float, default=None, help="Optional maximum Ecm filter.")
    ap.add_argument(
        "--plot",
        action="store_true",
        help="Also write a quick multi-panel plot of tracked/interpolated lines.",
    )
    ap.add_argument(
        "--plot-output",
        default=None,
        help="Optional plot file name. Default: <outdir>/qc_state_lines_continuous.png",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed parsing and tracking messages.",
    )
    return ap.parse_args()


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", s)


def collect_dat_files(inputs: List[str], workdir: Path, verbose: bool=False) -> List[Path]:
    files: List[Path] = []
    for raw in inputs:
        p = Path(raw)
        if not p.exists():
            print(f"[warn] input does not exist: {p}", file=sys.stderr)
            continue
        if p.is_file() and p.suffix.lower() == ".zip":
            dest = workdir / p.stem
            dest.mkdir(parents=True, exist_ok=True)
            if verbose:
                print(f"[read] extracting zip {p} -> {dest}")
            with zipfile.ZipFile(p, "r") as z:
                z.extractall(dest)
            for q in dest.rglob("*.dat"):
                files.append(q)
        elif p.is_dir():
            for q in p.rglob("*.dat"):
                files.append(q)
        elif p.is_file() and p.suffix.lower() == ".dat":
            files.append(p)
    # Filter: keep spectrum files, drop all_candidates and cache reports
    keep: List[Path] = []
    for f in files:
        name = f.name
        if "all_candidates" in name:
            continue
        if "QC_spectrum" in name or "spectrum_vs_L_all_levels" in name or "bestfit_QC_spectrum" in name:
            keep.append(f)
    keep = sorted(set(keep), key=lambda x: str(x))
    if verbose:
        print(f"[read] collected {len(keep)} spectrum .dat files")
        for f in keep:
            print("  ", f)
    return keep


def parse_spectrum_file(path: Path, irreps: set, lmin: float, lmax: float,
                        ymin: Optional[float], ymax: Optional[float]) -> List[Point]:
    pts: List[Point] = []
    with path.open("r", errors="replace") as f:
        for line_no, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 5:
                continue
            try:
                Lbyas = float(parts[0])
                spatial_L = float(parts[1])
                label = parts[2]
                state_index = int(float(parts[3]))
                Ecm = float(parts[4])
            except Exception:
                continue
            if label not in irreps:
                continue
            if Lbyas < lmin - 1e-9 or Lbyas > lmax + 1e-9:
                continue
            if ymin is not None and Ecm < ymin:
                continue
            if ymax is not None and Ecm > ymax:
                continue
            kind = "true_zero"
            if len(parts) >= 6:
                # v32x QC_spectrum has kind in col 5. v32y all_levels has no kind.
                if parts[5] in {"true_zero", "pole", "ZERO", "likely_zero"}:
                    kind = parts[5]
            # Keep true_zero-like rows only, if kind exists.
            if kind.lower() in {"pole", "likely_pole"}:
                continue
            pts.append(Point(Lbyas, spatial_L, label, Ecm, state_index, str(path), kind))
    return pts


def deduplicate_points(points: List[Point], tol: float=1e-10) -> List[Point]:
    # Deduplicate same label, same L, same Ecm from overlapping files.
    out: List[Point] = []
    seen = set()
    for p in sorted(points, key=lambda p: (p.label, p.Lbyas, p.Ecm, p.source)):
        key = (p.label, round(p.Lbyas, 9), round(p.Ecm / tol) if tol > 0 else p.Ecm)
        if key in seen:
            continue
        seen.add(key)
        out.append(p)
    return out


def group_by_label_L(points: List[Point]) -> Dict[str, Dict[float, List[Point]]]:
    d: Dict[str, Dict[float, List[Point]]] = {}
    for p in points:
        d.setdefault(p.label, {}).setdefault(p.Lbyas, []).append(p)
    for label in d:
        for L in d[label]:
            d[label][L].sort(key=lambda p: p.Ecm)
            # Assign rank after sorting; preserve original in state_index if needed.
            for i, p in enumerate(d[label][L]):
                p.state_index = i
    return d


def match_tracks_to_points(active: List[Track], new_points: List[Point], max_jump: float) -> Tuple[List[Tuple[int, int]], List[int], List[int]]:
    """Return matches as (track_index, point_index), unmatched track indices, unmatched point indices."""
    if not active or not new_points:
        return [], list(range(len(active))), list(range(len(new_points)))
    cost = np.full((len(active), len(new_points)), 1e9, dtype=float)
    for i, tr in enumerate(active):
        e0 = tr.last().Ecm
        for j, p in enumerate(new_points):
            c = abs(p.Ecm - e0)
            if c <= max_jump:
                cost[i, j] = c
    matches: List[Tuple[int, int]] = []
    if HAVE_SCIPY:
        rows, cols = linear_sum_assignment(cost)
        used_t, used_p = set(), set()
        for i, j in zip(rows, cols):
            if cost[i, j] <= max_jump:
                matches.append((int(i), int(j)))
                used_t.add(int(i)); used_p.add(int(j))
        unmatched_t = [i for i in range(len(active)) if i not in used_t]
        unmatched_p = [j for j in range(len(new_points)) if j not in used_p]
        return matches, unmatched_t, unmatched_p
    # Fallback greedy by sorted edges
    edges = []
    for i in range(len(active)):
        for j in range(len(new_points)):
            if cost[i, j] <= max_jump:
                edges.append((cost[i, j], i, j))
    edges.sort()
    used_t, used_p = set(), set()
    for _, i, j in edges:
        if i not in used_t and j not in used_p:
            matches.append((i, j))
            used_t.add(i); used_p.add(j)
    unmatched_t = [i for i in range(len(active)) if i not in used_t]
    unmatched_p = [j for j in range(len(new_points)) if j not in used_p]
    return matches, unmatched_t, unmatched_p


def track_label(label: str, by_L: Dict[float, List[Point]], max_jump: float, max_gap: float) -> List[Track]:
    Ls = sorted(by_L.keys())
    tracks: List[Track] = []
    active: List[Track] = []
    next_id = 0
    last_L: Optional[float] = None
    for L in Ls:
        pts = by_L[L]
        if last_L is None or abs(L - last_L) > max_gap:
            # Close active tracks and start anew after a big gap.
            tracks.extend(active)
            active = []
        matches, unmatched_t, unmatched_p = match_tracks_to_points(active, pts, max_jump=max_jump)
        new_active: List[Track] = []
        matched_t = set()
        matched_p = set()
        for ti, pj in matches:
            tr = active[ti]
            tr.add(pts[pj])
            new_active.append(tr)
            matched_t.add(ti); matched_p.add(pj)
        # Close unmatched old active tracks.
        for ti, tr in enumerate(active):
            if ti not in matched_t:
                tracks.append(tr)
        # New tracks from unmatched points.
        for pj, p in enumerate(pts):
            if pj not in matched_p:
                tr = Track(label=label, track_id=next_id, points=[p])
                next_id += 1
                new_active.append(tr)
        # Sort active by last energy for determinism.
        active = sorted(new_active, key=lambda tr: tr.last().Ecm)
        last_L = L
    tracks.extend(active)
    # Renumber by first L then first E
    tracks = sorted(tracks, key=lambda tr: (tr.sorted_points()[0].Lbyas, tr.sorted_points()[0].Ecm, tr.track_id))
    for i, tr in enumerate(tracks):
        tr.track_id = i
    return tracks


def interp_track(track: Track, method: str, dense_n: int, xi: float) -> Optional[np.ndarray]:
    pts = track.sorted_points()
    # average duplicate L, if any
    L_to_E: Dict[float, List[float]] = {}
    for p in pts:
        L_to_E.setdefault(float(p.Lbyas), []).append(float(p.Ecm))
    xs = np.array(sorted(L_to_E.keys()), dtype=float)
    ys = np.array([np.mean(L_to_E[x]) for x in xs], dtype=float)
    if len(xs) < 2:
        return None
    L_dense = np.linspace(xs.min(), xs.max(), dense_n)
    actual_method = method
    if not HAVE_SCIPY:
        actual_method = "linear"
        E_dense = np.interp(L_dense, xs, ys)
    elif method == "linear" or len(xs) < 4:
        # cubic requires enough points; use linear for short tracks.
        actual_method = "linear" if method != "linear" else method
        f = interp1d(xs, ys, kind="linear", bounds_error=False, fill_value=np.nan)
        E_dense = f(L_dense)
    elif method == "cubic":
        f = CubicSpline(xs, ys, bc_type="natural")
        E_dense = f(L_dense)
    elif method == "pchip":
        f = PchipInterpolator(xs, ys)
        E_dense = f(L_dense)
    else:
        E_dense = np.interp(L_dense, xs, ys)
        actual_method = "linear"
    spatial = L_dense * xi
    return np.column_stack([L_dense, spatial, E_dense, np.full_like(L_dense, len(xs), dtype=float)])


def write_outputs(tracks_by_label: Dict[str, List[Track]], outdir: Path, args: argparse.Namespace) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    raw_dir = outdir / "tracked_points"
    dense_dir = outdir / f"interpolated_{args.interp}"
    raw_dir.mkdir(parents=True, exist_ok=True)
    dense_dir.mkdir(parents=True, exist_ok=True)

    all_raw_lines = ["# label track_id npoints Lbyas spatial_L Ecm source_state_index source_file\n"]
    all_dense_lines = ["# label track_id interp_method n_observed Lbyas spatial_L Ecm\n"]
    summary_lines = [
        "# label track_id npoints Lmin_obs Lmax_obs Emin_obs Emax_obs wrote_interpolated\n"
    ]

    for label, tracks in tracks_by_label.items():
        for tr in tracks:
            pts = tr.sorted_points()
            npts = len({p.Lbyas for p in pts})
            if npts < args.min_points:
                wrote = 0
                summary_lines.append(
                    f"{label} {tr.track_id} {npts} {pts[0].Lbyas:.12g} {pts[-1].Lbyas:.12g} {min(p.Ecm for p in pts):.17g} {max(p.Ecm for p in pts):.17g} {wrote}\n"
                )
                continue
            if args.require_full_span:
                Lset = {round(p.Lbyas, 9) for p in pts}
                if round(args.lmin, 9) not in Lset or round(args.lmax, 9) not in Lset:
                    wrote = 0
                    summary_lines.append(
                        f"{label} {tr.track_id} {npts} {pts[0].Lbyas:.12g} {pts[-1].Lbyas:.12g} {min(p.Ecm for p in pts):.17g} {max(p.Ecm for p in pts):.17g} {wrote}\n"
                    )
                    continue
            raw_name = raw_dir / f"{safe_name(label)}_state_{tr.track_id:03d}_tracked_points.dat"
            with raw_name.open("w") as f:
                f.write("# Lbyas spatial_L Ecm label track_id source_state_index source_file\n")
                for p in pts:
                    f.write(f"{p.Lbyas:.12g} {p.spatial_L:.17g} {p.Ecm:.17g} {label} {tr.track_id} {p.state_index} {p.source}\n")
                    all_raw_lines.append(
                        f"{label} {tr.track_id} {npts} {p.Lbyas:.12g} {p.spatial_L:.17g} {p.Ecm:.17g} {p.state_index} {p.source}\n"
                    )
            arr = interp_track(tr, args.interp, args.dense_n, args.xi)
            if arr is None:
                wrote = 0
            else:
                wrote = 1
                dense_name = dense_dir / f"{safe_name(label)}_state_{tr.track_id:03d}_interp_{args.interp}.dat"
                # actual method could be linear for short cubic/pchip tracks; infer by points.
                actual = args.interp
                if args.interp in {"cubic", "pchip"} and len({p.Lbyas for p in pts}) < 4:
                    actual = "linear"
                with dense_name.open("w") as f:
                    f.write("# Lbyas spatial_L Ecm label track_id interp_method n_observed\n")
                    for Lbyas, spatial, Ecm, nobsv in arr:
                        f.write(f"{Lbyas:.12g} {spatial:.17g} {Ecm:.17g} {label} {tr.track_id} {actual} {int(nobsv)}\n")
                        all_dense_lines.append(
                            f"{label} {tr.track_id} {actual} {int(nobsv)} {Lbyas:.12g} {spatial:.17g} {Ecm:.17g}\n"
                        )
            summary_lines.append(
                f"{label} {tr.track_id} {npts} {pts[0].Lbyas:.12g} {pts[-1].Lbyas:.12g} {min(p.Ecm for p in pts):.17g} {max(p.Ecm for p in pts):.17g} {wrote}\n"
            )

    (outdir / "all_tracked_points.dat").write_text("".join(all_raw_lines))
    (outdir / f"all_interpolated_state_lines_{args.interp}.dat").write_text("".join(all_dense_lines))
    (outdir / "state_line_summary.dat").write_text("".join(summary_lines))


def make_plot(outdir: Path, args: argparse.Namespace, tracks_by_label: Dict[str, List[Track]]) -> None:
    import matplotlib.pyplot as plt
    irrep_list = [r for r in args.irreps if r in tracks_by_label]
    if not irrep_list:
        return
    fig, axes = plt.subplots(1, len(irrep_list), figsize=(4.2 * len(irrep_list), 7), sharey=True)
    if len(irrep_list) == 1:
        axes = [axes]
    dense_file = outdir / f"all_interpolated_state_lines_{args.interp}.dat"
    data_by = {}
    with dense_file.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split()
            label = parts[0]; track = int(parts[1]); L = float(parts[4]); E = float(parts[6])
            data_by.setdefault((label, track), []).append((L, E))
    raw_by = {}
    with (outdir / "all_tracked_points.dat").open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split()
            label = parts[0]; track = int(parts[1]); L = float(parts[3]); E = float(parts[5])
            raw_by.setdefault((label, track), []).append((L, E))
    kkpi = 0.26302
    kkpipi = 0.33208
    for ax, label in zip(axes, irrep_list):
        for (lab, track), vals in data_by.items():
            if lab != label:
                continue
            arr = np.array(vals)
            ax.plot(arr[:,0], arr[:,1], '-', lw=1.2)
        for (lab, track), vals in raw_by.items():
            if lab != label:
                continue
            arr = np.array(vals)
            ax.plot(arr[:,0], arr[:,1], 'o', ms=3)
        ax.set_title(label)
        ax.set_xlim(args.lmin, args.lmax)
        ax.axhline(kkpi, color='0.7', ls='--', lw=0.8)
        ax.axhline(kkpipi, color='0.7', ls='--', lw=0.8)
        ax.set_xlabel(r'$L/a_s$')
        if args.ymin is not None or args.ymax is not None:
            ax.set_ylim(args.ymin, args.ymax)
    axes[0].set_ylabel(r'$a_t E_{\rm cm}$')
    fig.tight_layout()
    out = Path(args.plot_output) if args.plot_output else outdir / "qc_state_lines_continuous.png"
    fig.savefig(out, dpi=200)
    print(f"[plot] wrote {out}")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    tmp = Path(tempfile.mkdtemp(prefix="qc_state_lines_extract_"))
    try:
        files = collect_dat_files(args.inputs, tmp, verbose=args.verbose)
        if not files:
            raise SystemExit("[error] no spectrum .dat files found in inputs")
        irreps = set(args.irreps)
        points: List[Point] = []
        for f in files:
            pts = parse_spectrum_file(f, irreps, args.lmin, args.lmax, args.ymin, args.ymax)
            if args.verbose:
                print(f"[parse] {f}: {len(pts)} points")
            points.extend(pts)
        points = deduplicate_points(points)
        if not points:
            raise SystemExit("[error] no QC spectrum points found after filters")
        print(f"[info] parsed {len(points)} unique QC points")
        by_label_L = group_by_label_L(points)
        tracks_by_label: Dict[str, List[Track]] = {}
        for label in args.irreps:
            if label not in by_label_L:
                print(f"[warn] no points for {label}")
                continue
            tracks = track_label(label, by_label_L[label], max_jump=args.max_jump, max_gap=args.max_gap)
            tracks_by_label[label] = tracks
            n_written = sum(1 for tr in tracks if len({p.Lbyas for p in tr.points}) >= args.min_points)
            print(f"[track] {label}: {len(tracks)} raw tracks, {n_written} with >= {args.min_points} points")
        write_outputs(tracks_by_label, outdir, args)
        print(f"[write] {outdir / 'state_line_summary.dat'}")
        print(f"[write] {outdir / 'all_tracked_points.dat'}")
        print(f"[write] {outdir / ('all_interpolated_state_lines_' + args.interp + '.dat')}")
        print(f"[write] per-state files under {outdir / ('interpolated_' + args.interp)}")
        if args.plot:
            make_plot(outdir, args, tracks_by_label)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
