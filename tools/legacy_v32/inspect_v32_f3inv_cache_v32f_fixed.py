#!/usr/bin/env python3
"""
Robust inspector for v32 F3inv/Vsel binary cache files.

It reads the v32f-style cache metadata (success is a 4-byte int, matching ProjectedQCCacheEntry::success) without loading the large matrices and reports:
  - complete entries per irrep
  - unique Ecm values per irrep
  - coarse-grid completeness for a requested Ecm range/coarseN
  - extra/refined Ecm points
  - whether the binary file appears truncated/corrupted

Supported magic:
  V32F_F3INV_FULL_CACHE_1
"""
import argparse
import os
import struct
import sys
from bisect import bisect_left
from collections import defaultdict

U64 = "<Q"
I32 = "<i"
F64 = "<d"
BOOL = "<?"

class CacheEOF(Exception):
    pass

class CacheFormatError(Exception):
    pass


def file_size(path):
    return os.path.getsize(path)


def read_exact(f, n, what="data"):
    pos = f.tell()
    b = f.read(n)
    if len(b) != n:
        raise CacheEOF(f"Unexpected EOF while reading {what}: offset={pos} wanted={n} got={len(b)}")
    return b


def read_u64(f):
    return struct.unpack(U64, read_exact(f, 8, "uint64"))[0]


def read_i32(f):
    return struct.unpack(I32, read_exact(f, 4, "int32"))[0]


def read_f64(f):
    return struct.unpack(F64, read_exact(f, 8, "double"))[0]


def read_bool(f):
    return struct.unpack(BOOL, read_exact(f, 1, "bool"))[0]


def read_string(f):
    n = read_u64(f)
    if n > 100_000_000:
        raise CacheFormatError(f"Unreasonable string length {n} at offset={f.tell()-8}; wrong format or corrupted cache")
    if n == 0:
        return ""
    return read_exact(f, n, "string").decode("utf-8", errors="replace")


def seek_checked(f, nbytes, what, fsize):
    pos = f.tell()
    if nbytes < 0:
        raise CacheFormatError(f"Negative skip {nbytes} for {what} at offset={pos}")
    if pos + nbytes > fsize:
        raise CacheEOF(f"Unexpected EOF while skipping {what}: offset={pos} wanted_skip={nbytes} file_size={fsize}")
    f.seek(nbytes, os.SEEK_CUR)


def skip_comp(f, fsize):
    seek_checked(f, 16, "complex<double>", fsize)


def skip_vec_comp(f, fsize):
    n = read_u64(f)
    seek_checked(f, 16 * n, f"vector<complex<double>> length={n}", fsize)


def skip_nested_vec_comp(f, fsize):
    n = read_u64(f)
    for _ in range(n):
        skip_vec_comp(f, fsize)


def skip_matrix(f, fsize):
    rows = read_u64(f)
    cols = read_u64(f)
    if rows > 2_000_000 or cols > 2_000_000:
        raise CacheFormatError(f"Unreasonable matrix shape ({rows},{cols}) at offset={f.tell()-16}; wrong format/corrupted")
    nbytes = 16 * rows * cols
    seek_checked(f, nbytes, f"MatrixXcd shape=({rows},{cols})", fsize)
    return rows, cols


def parse_cache(path, only_irrep=None, keep_errors=False, verbose=False):
    entries = []
    problems = []
    fsize = file_size(path)
    with open(path, "rb") as f:
        try:
            magic = read_string(f)
            if magic != "V32F_F3INV_FULL_CACHE_1":
                raise CacheFormatError(f"Unsupported magic/version: {magic!r}")
            ncache = read_u64(f)
            if verbose:
                print(f"[cache] magic={magic} ncache={ncache} size={fsize}")

            for ci in range(ncache):
                block_start = f.tell()
                try:
                    cache_label = read_string(f)
                    ngrid = read_u64(f)
                except (CacheEOF, CacheFormatError) as e:
                    problems.append(f"block {ci}: failed header at offset={block_start}: {e}")
                    break

                if verbose:
                    print(f"[cache] block={ci} label={cache_label} expected_ngrid={ngrid} offset={block_start}")

                complete_in_block = 0
                for gi in range(ngrid):
                    entry_start = f.tell()
                    try:
                        label = read_string(f)
                        idx = read_i32(f)
                        Ecm = read_f64(f)
                        En = read_f64(f)
                        success = read_i32(f)
                        total_dim = read_i32(f)
                        proj_dim = read_i32(f)

                        F3proj_shape = skip_matrix(f, fsize)
                        F3full_shape = skip_matrix(f, fsize)
                        skip_nested_vec_comp(f, fsize)  # plm_config
                        skip_nested_vec_comp(f, fsize)  # klm_config
                        skip_vec_comp(f, fsize)         # total_P
                        skip_comp(f, fsize)             # En_c
                        V_shape = skip_matrix(f, fsize)
                        projector_idem = read_f64(f)
                        projector_closure = read_f64(f)
                        F3_equiv = read_f64(f)
                        F3_leak = read_f64(f)
                        error = read_string(f)
                        entry_end = f.tell()
                    except (CacheEOF, CacheFormatError) as e:
                        problems.append(
                            f"block {ci} label={cache_label} entry {gi}/{ngrid}: incomplete/corrupt at offset={entry_start}: {e}"
                        )
                        # Stop parsing: any later bytes are not reliable because we do not know alignment.
                        return entries, problems, fsize

                    complete_in_block += 1
                    if only_irrep is None or label == only_irrep or cache_label == only_irrep:
                        if keep_errors or success or True:
                            entries.append({
                                "cache_label": cache_label,
                                "label": label,
                                "i": idx,
                                "Ecm": Ecm,
                                "En": En,
                                "success": bool(success),
                                "total_dim": total_dim,
                                "proj_dim": proj_dim,
                                "F3inv_proj_shape": F3proj_shape,
                                "F3inv_full_shape": F3full_shape,
                                "Vsel_shape": V_shape,
                                "projector_idem": projector_idem,
                                "projector_closure": projector_closure,
                                "F3_equiv": F3_equiv,
                                "F3_leak": F3_leak,
                                "error": error,
                                "byte_start": entry_start,
                                "byte_end": entry_end,
                            })

                if verbose:
                    print(f"[cache] block={ci} label={cache_label} complete_entries={complete_in_block}/{ngrid}")
        except (CacheEOF, CacheFormatError) as e:
            problems.append(f"file header/top-level parse failed: {e}")
    return entries, problems, fsize


def unique_sorted(values, tol=1e-12):
    vals = sorted(float(x) for x in values if x == x)
    if not vals:
        return []
    out = [vals[0]]
    for x in vals[1:]:
        if abs(x - out[-1]) > tol:
            out.append(x)
    return out


def has_near(sorted_vals, x, tol):
    j = bisect_left(sorted_vals, x)
    cand = []
    if j < len(sorted_vals):
        cand.append(sorted_vals[j])
    if j > 0:
        cand.append(sorted_vals[j-1])
    return any(abs(y - x) <= tol for y in cand)


def expected_coarse_grid(coarse_n, ecm_min, ecm_max):
    if coarse_n <= 1:
        return []
    return [ecm_min + (ecm_max - ecm_min) * k / (coarse_n - 1) for k in range(coarse_n)]


def summarize(entries, coarse_n=None, ecm_min=None, ecm_max=None, tol=1e-10, max_print=20):
    groups = defaultdict(list)
    for e in entries:
        groups[e["label"]].append(e)
    print("\n========== cache summary ==========")
    print(f"complete entries parsed = {len(entries)}")
    if not groups:
        print("No complete entries parsed.")
        return

    for label in sorted(groups):
        es = groups[label]
        energies = unique_sorted([x["Ecm"] for x in es], tol=tol/10)
        nsuccess = sum(1 for x in es if x["success"])
        print(f"\n[{label}]")
        print(f"  entries          = {len(es)}")
        print(f"  unique_Ecm       = {len(energies)}")
        print(f"  success_entries  = {nsuccess}")
        print(f"  failed_entries   = {len(es)-nsuccess}")
        print(f"  Ecm_min/max      = {energies[0]:.17g}  {energies[-1]:.17g}" if energies else "  Ecm_min/max      = none")
        shapes = defaultdict(int)
        for e in es:
            shapes[(e["F3inv_full_shape"], e["Vsel_shape"])] += 1
        for (f3shape, vshape), cnt in sorted(shapes.items(), key=lambda kv: -kv[1])[:5]:
            print(f"  shape count      = {cnt:8d}  F3inv_full={f3shape} Vsel={vshape}")

        if energies:
            show = energies[:max_print]
            print("  first energies   = " + " ".join(f"{x:.12f}" for x in show))
            if len(energies) > max_print:
                print("  last energies    = " + " ".join(f"{x:.12f}" for x in energies[-max_print:]))

        if coarse_n and ecm_min is not None and ecm_max is not None:
            expected = expected_coarse_grid(coarse_n, ecm_min, ecm_max)
            present = [x for x in expected if has_near(energies, x, tol)]
            missing = [x for x in expected if not has_near(energies, x, tol)]
            extra = [x for x in energies if not has_near(expected, x, tol)]
            print(f"  expected_coarseN = {coarse_n}")
            print(f"  coarse_present   = {len(present)}")
            print(f"  missing_coarse   = {len(missing)}")
            print(f"  refined_or_extra = {len(extra)}")
            if missing:
                print("  first missing    = " + " ".join(f"{x:.12f}" for x in missing[:max_print]))
            if extra:
                print("  first extra      = " + " ".join(f"{x:.12f}" for x in extra[:max_print]))


def write_energies(path, entries, tol=1e-12):
    groups = defaultdict(list)
    for e in entries:
        groups[e["label"]].append(e["Ecm"])
    with open(path, "w") as f:
        f.write("# label Ecm\n")
        for label in sorted(groups):
            for E in unique_sorted(groups[label], tol=tol):
                f.write(f"{label} {E:.17g}\n")
    print(f"[write] energies written to {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True, help="binary cache file")
    ap.add_argument("--irrep", default=None, help="optional irrep label, e.g. 110_A2")
    ap.add_argument("--coarse-n", type=int, default=None)
    ap.add_argument("--ecm-min", type=float, default=None)
    ap.add_argument("--ecm-max", type=float, default=None)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--max-print", type=int, default=20)
    ap.add_argument("--write-energies", default=None)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(args.cache):
        print(f"[error] cache file does not exist: {args.cache}", file=sys.stderr)
        sys.exit(2)

    print(f"[inspect] cache = {args.cache}")
    print(f"[inspect] size  = {file_size(args.cache)} bytes")
    entries, problems, fsize = parse_cache(args.cache, only_irrep=args.irrep, verbose=args.verbose)
    if problems:
        print("\n========== parse warnings/problems ==========")
        for p in problems:
            print("[problem] " + p)
        print("[problem] The file may be truncated/incomplete, or this script's format does not match the cache exactly.")
        print("[problem] Summary below uses only complete entries parsed before the failure.\n")

    summarize(entries, coarse_n=args.coarse_n, ecm_min=args.ecm_min, ecm_max=args.ecm_max, tol=args.tol, max_print=args.max_print)

    if args.write_energies:
        write_energies(args.write_energies, entries, tol=args.tol/10)

    if problems:
        sys.exit(1)

if __name__ == "__main__":
    main()
