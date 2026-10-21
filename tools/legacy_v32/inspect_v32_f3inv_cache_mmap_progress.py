#!/usr/bin/env python3
"""
Fast mmap/progress inspector for v32x/v32f F3inv/Vsel binary cache files.

Purpose:
  - Print how many Ecm values are stored per irrep/block.
  - Print/write the stored energies.
  - Check whether the initial coarse grid of coarseN points is complete.
  - Count refined/extra energies beyond the expected coarse grid.

This parser matches the v32f/v32x binary format where success is stored as a
4-byte int, not a 1-byte bool.

Expected magic:
  V32F_F3INV_FULL_CACHE_1

It does NOT load matrix values; it only skips over them using the stored matrix
row/column dimensions. It uses mmap so it is much faster than repeated reads on
multi-GB cache files.
"""
import argparse
import mmap
import os
import struct
import sys
from bisect import bisect_left
from collections import defaultdict

U64 = struct.Struct('<Q')
I32 = struct.Struct('<i')
F64 = struct.Struct('<d')

class ParseError(Exception):
    pass


def need(mm, off, n, what):
    if off + n > len(mm):
        raise ParseError(f"EOF while reading/skipping {what}: offset={off} need={n} file_size={len(mm)}")


def u64(mm, off):
    need(mm, off, 8, 'uint64')
    return U64.unpack_from(mm, off)[0], off + 8


def i32(mm, off):
    need(mm, off, 4, 'int32')
    return I32.unpack_from(mm, off)[0], off + 4


def f64(mm, off):
    need(mm, off, 8, 'double')
    return F64.unpack_from(mm, off)[0], off + 8


def string(mm, off, decode=True):
    n, off = u64(mm, off)
    if n > 100_000_000:
        raise ParseError(f"unreasonable string length {n} at offset={off-8}")
    need(mm, off, n, f'string length={n}')
    if decode:
        s = mm[off:off+n].decode('utf-8', errors='replace')
    else:
        s = None
    return s, off + n


def skip_comp(mm, off):
    need(mm, off, 16, 'complex<double>')
    return off + 16


def skip_vec_comp(mm, off):
    n, off = u64(mm, off)
    nbytes = 16 * n
    need(mm, off, nbytes, f'vector<complex<double>> length={n}')
    return off + nbytes


def skip_nested_vec_comp(mm, off):
    nouter, off = u64(mm, off)
    # The momentum/config arrays are usually small-ish. We still must read each
    # inner length because the vectors are jagged.
    for _ in range(nouter):
        off = skip_vec_comp(mm, off)
    return off


def skip_matrix(mm, off, max_dim=2_000_000):
    rows, off = u64(mm, off)
    cols, off = u64(mm, off)
    if rows > max_dim or cols > max_dim:
        raise ParseError(f"unreasonable matrix shape ({rows},{cols}) at offset={off-16}; wrong format/alignment")
    nbytes = 16 * rows * cols
    need(mm, off, nbytes, f'MatrixXcd shape=({rows},{cols})')
    return (rows, cols), off + nbytes


def expected_grid(coarse_n, e_min, e_max):
    if coarse_n is None or coarse_n <= 1 or e_min is None or e_max is None:
        return []
    return [e_min + (e_max - e_min) * k / (coarse_n - 1) for k in range(coarse_n)]


def unique_sorted(vals, tol=1e-12):
    vals = sorted(float(v) for v in vals)
    if not vals:
        return []
    out = [vals[0]]
    for v in vals[1:]:
        if abs(v - out[-1]) > tol:
            out.append(v)
    return out


def has_near(sorted_vals, x, tol):
    j = bisect_left(sorted_vals, x)
    if j < len(sorted_vals) and abs(sorted_vals[j] - x) <= tol:
        return True
    if j > 0 and abs(sorted_vals[j-1] - x) <= tol:
        return True
    return False


def parse_one_entry(mm, off, keep_entry=False):
    start = off
    label, off = string(mm, off, decode=True)
    idx, off = i32(mm, off)
    Ecm, off = f64(mm, off)
    En, off = f64(mm, off)
    success, off = i32(mm, off)  # IMPORTANT: v32 cache writes int success, 4 bytes
    total_dim, off = i32(mm, off)
    proj_dim, off = i32(mm, off)

    F3proj_shape, off = skip_matrix(mm, off)
    F3full_shape, off = skip_matrix(mm, off)
    off = skip_nested_vec_comp(mm, off)  # plm_config
    off = skip_nested_vec_comp(mm, off)  # klm_config
    off = skip_vec_comp(mm, off)         # total_P
    off = skip_comp(mm, off)             # En_c
    V_shape, off = skip_matrix(mm, off)
    projector_idem, off = f64(mm, off)
    projector_closure, off = f64(mm, off)
    F3_equiv, off = f64(mm, off)
    F3_leak, off = f64(mm, off)
    # Decode error only if keeping details; otherwise just skip string.
    error, off = string(mm, off, decode=keep_entry)

    if keep_entry:
        return {
            'label': label,
            'idx': idx,
            'Ecm': Ecm,
            'En': En,
            'success': bool(success),
            'total_dim': total_dim,
            'proj_dim': proj_dim,
            'F3inv_proj_shape': F3proj_shape,
            'F3inv_full_shape': F3full_shape,
            'Vsel_shape': V_shape,
            'projector_idem': projector_idem,
            'projector_closure': projector_closure,
            'F3_equiv': F3_equiv,
            'F3_leak': F3_leak,
            'error': error if error is not None else '',
            'byte_start': start,
            'byte_end': off,
        }, off
    else:
        # Minimal tuple for speed/memory.
        return (label, idx, Ecm, En, bool(success), total_dim, proj_dim, start, off), off


def inspect_cache(path, only_irrep=None, progress_every=500, max_entries=None, quiet=False):
    size = os.path.getsize(path)
    if not quiet:
        print(f"[inspect] cache = {path}", flush=True)
        print(f"[inspect] size  = {size} bytes", flush=True)
        print("[inspect] mmap opening...", flush=True)

    summaries = defaultdict(lambda: {
        'entries': 0,
        'success': 0,
        'failed': 0,
        'energies': [],
        'dims': defaultdict(int),
        'byte_min': None,
        'byte_max': None,
    })
    blocks = []
    problems = []

    with open(path, 'rb') as fh:
        with mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            off = 0
            magic, off = string(mm, off, decode=True)
            if magic != 'V32F_F3INV_FULL_CACHE_1':
                raise ParseError(f"unsupported magic {magic!r}; expected V32F_F3INV_FULL_CACHE_1")
            ncache, off = u64(mm, off)
            if not quiet:
                print(f"[inspect] magic={magic} ncache_blocks={ncache}", flush=True)

            total_seen = 0
            for bi in range(ncache):
                block_off = off
                cache_label, off = string(mm, off, decode=True)
                ngrid, off = u64(mm, off)
                blocks.append((bi, cache_label, ngrid, block_off))
                if not quiet:
                    print(f"[block {bi}] label={cache_label} entries={ngrid} offset={block_off}", flush=True)

                for gi in range(ngrid):
                    entry_off = off
                    try:
                        item, off = parse_one_entry(mm, off, keep_entry=False)
                    except Exception as e:
                        problems.append(f"block={bi} label={cache_label} entry={gi}/{ngrid} offset={entry_off}: {e}")
                        return summaries, blocks, problems, size, off

                    label, idx, Ecm, En, success, total_dim, proj_dim, start, end = item
                    # Normally label == cache_label. Use either, but summarize by actual entry label.
                    if only_irrep is None or label == only_irrep or cache_label == only_irrep:
                        s = summaries[label]
                        s['entries'] += 1
                        s['success'] += int(success)
                        s['failed'] += int(not success)
                        s['energies'].append(Ecm)
                        s['dims'][(total_dim, proj_dim)] += 1
                        s['byte_min'] = start if s['byte_min'] is None else min(s['byte_min'], start)
                        s['byte_max'] = end if s['byte_max'] is None else max(s['byte_max'], end)

                    total_seen += 1
                    if (not quiet) and progress_every and ((gi + 1) % progress_every == 0 or gi + 1 == ngrid):
                        pct = 100.0 * (gi + 1) / max(1, ngrid)
                        print(f"[block {bi} {cache_label}] {gi+1}/{ngrid} ({pct:.1f}%) offset={off}", flush=True)

                    if max_entries is not None and total_seen >= max_entries:
                        if not quiet:
                            print(f"[inspect] stopped early at max_entries={max_entries}", flush=True)
                        return summaries, blocks, problems, size, off

            return summaries, blocks, problems, size, off


def print_summary(summaries, coarse_n, ecm_min, ecm_max, tol, max_print, write_energies):
    print("\n========== cache summary ==========")
    if not summaries:
        print("No matching entries summarized.")
        return

    expected = expected_grid(coarse_n, ecm_min, ecm_max)

    write_lines = []
    for label in sorted(summaries.keys()):
        s = summaries[label]
        uniq = unique_sorted(s['energies'], tol=1e-12)
        print(f"\n[{label}]")
        print(f"entries              = {s['entries']}")
        print(f"success_entries      = {s['success']}")
        print(f"failed_entries       = {s['failed']}")
        print(f"unique_Ecm           = {len(uniq)}")
        if uniq:
            print(f"Ecm_min              = {uniq[0]:.15f}")
            print(f"Ecm_max              = {uniq[-1]:.15f}")
        print(f"byte_range           = {s['byte_min']} .. {s['byte_max']}")
        if s['dims']:
            dim_txt = ', '.join([f"total={k[0]},proj={k[1]}: {v}" for k, v in sorted(s['dims'].items())])
            print(f"dims(counts)         = {dim_txt}")

        if expected:
            n_present = sum(1 for x in expected if has_near(uniq, x, tol))
            n_missing = len(expected) - n_present
            n_extra = sum(1 for x in uniq if not has_near(expected, x, tol))
            print(f"coarse_points_present= {n_present}/{len(expected)}")
            print(f"missing_coarse_points= {n_missing}")
            print(f"refined_or_extra     = {n_extra}")
            if n_missing:
                missing = [x for x in expected if not has_near(uniq, x, tol)]
                print("first_missing_Ecm    = " + ' '.join(f"{x:.15f}" for x in missing[:max_print]))
            if n_extra:
                extra = [x for x in uniq if not has_near(expected, x, tol)]
                print("first_extra_Ecm      = " + ' '.join(f"{x:.15f}" for x in extra[:max_print]))

        print("first_Ecm_values     = " + ' '.join(f"{x:.15f}" for x in uniq[:max_print]))
        print("last_Ecm_values      = " + ' '.join(f"{x:.15f}" for x in uniq[-max_print:]))

        for E in uniq:
            write_lines.append(f"{label} {E:.17g}\n")

    if write_energies:
        with open(write_energies, 'w') as f:
            f.write("# label Ecm\n")
            f.writelines(write_lines)
        print(f"\n[write] energies -> {write_energies}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache', required=True, help='Path to *_F3inv_Vsel_merged.bin cache file')
    ap.add_argument('--irrep', default=None, help='Only summarize this irrep, e.g. 110_A2')
    ap.add_argument('--coarse-n', type=int, default=None, help='Expected coarseN, e.g. 10000')
    ap.add_argument('--ecm-min', type=float, default=None, help='Expected coarse Ecm minimum')
    ap.add_argument('--ecm-max', type=float, default=None, help='Expected coarse Ecm maximum')
    ap.add_argument('--tol', type=float, default=1e-10, help='Tolerance for matching expected coarse-grid energies')
    ap.add_argument('--write-energies', default=None, help='Write unique Ecm values to this text file')
    ap.add_argument('--progress-every', type=int, default=250, help='Print progress every N entries per block; 0 disables')
    ap.add_argument('--max-print', type=int, default=20, help='How many first/last/missing/extra energies to print')
    ap.add_argument('--max-entries', type=int, default=None, help='Debug: stop after parsing this many entries total')
    ap.add_argument('--quiet', action='store_true', help='Suppress progress prints')
    args = ap.parse_args()

    try:
        summaries, blocks, problems, size, final_off = inspect_cache(
            args.cache,
            only_irrep=args.irrep,
            progress_every=args.progress_every,
            max_entries=args.max_entries,
            quiet=args.quiet,
        )
    except Exception as e:
        print(f"[fatal] {e}", file=sys.stderr)
        sys.exit(2)

    if problems:
        print("\n========== parse warnings/problems ==========")
        for p in problems:
            print("[problem] " + p)
        print("[problem] Summary below uses only entries parsed before the failure.")

    print(f"\n[inspect] final_offset={final_off} file_size={size} parsed_fraction={final_off/max(1,size):.6f}")
    print_summary(summaries, args.coarse_n, args.ecm_min, args.ecm_max, args.tol, args.max_print, args.write_energies)


if __name__ == '__main__':
    main()
