#!/usr/bin/env python3
"""
Compare two CFE discharge time series and report statistics.

Handles three output formats:
  - ngen CSV:       header row, columns: step, time, discharge, ...
  - standalone:     # comment lines, then: timestep discharge
  - fred CSV:       2 header lines, columns: datetime, discharge

Usage:
    compare_discharge.py <reference> <test_output> [--tol 1e-6] [--ref-fmt auto] [--test-fmt auto]

Exit code: 0 = PASS, 1 = FAIL, 2 = input error
"""

import argparse
import csv
import os
import sys


def read_ngen_csv(path):
    """ngen output: 1 header, columns: step, time, value, ..."""
    values = {}
    with open(path) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            values[int(row[0])] = float(row[2])
    return values


def read_standalone(path):
    """cfe_bmi_driver output: # comments, then space-delimited timestep value."""
    values = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                values[int(parts[0])] = float(parts[1])
    return values


def read_fred_csv(path):
    """Fred's output: 2 header lines, columns: datetime, discharge."""
    values = {}
    with open(path) as f:
        reader = csv.reader(f)
        next(reader)
        next(reader)
        for i, row in enumerate(reader):
            values[i] = float(row[1])
    return values


def detect_format(path):
    """Guess file format from first few lines."""
    with open(path) as f:
        lines = []
        for line in f:
            lines.append(line.strip())
            if len(lines) >= 5:
                break

    if not lines:
        return "ngen"

    # Find the first non-comment, non-empty line
    data_line = None
    comment_count = 0
    for line in lines:
        if line.startswith("#"):
            comment_count += 1
            continue
        if line:
            data_line = line
            break

    if data_line is None:
        return "ngen"

    # Files starting with # comments: distinguish fred CSV from standalone
    if comment_count > 0:
        if "," in data_line:
            return "fred"
        return "standalone"

    # No comment lines: ngen CSV has a text header, data on line 2
    parts = lines[0].split(",")
    if len(parts) >= 3:
        if len(lines) > 1:
            second_parts = lines[1].split(",")
            if second_parts[0].strip().isdigit():
                return "ngen"
        return "fred"

    return "ngen"


READERS = {
    "ngen": read_ngen_csv,
    "standalone": read_standalone,
    "fred": read_fred_csv,
}


def compare(ref_q, test_q, abs_threshold, trunc_bound=5e-9):
    """Compare two discharge dicts {step: value}. Returns (passed, stats)."""
    steps = sorted(set(ref_q.keys()) & set(test_q.keys()))
    n = len(steps)
    if n == 0:
        return False, {"error": "no overlapping timesteps"}

    exact = 0
    within_trunc = 0
    beyond_trunc = 0
    max_abs = 0.0
    max_rel = 0.0
    max_abs_step = 0
    max_rel_step = 0

    for step in steps:
        rq = ref_q[step]
        tq = test_q[step]
        ad = abs(rq - tq)

        if ad == 0:
            exact += 1
        elif ad <= trunc_bound:
            within_trunc += 1
        else:
            beyond_trunc += 1

        if ad > max_abs:
            max_abs = ad
            max_abs_step = step
        if rq != 0:
            rd = ad / abs(rq) * 100
            if rd > max_rel:
                max_rel = rd
                max_rel_step = step

    total = exact + within_trunc + beyond_trunc
    passed = max_abs < abs_threshold

    stats = {
        "total": total,
        "exact": exact,
        "within_trunc": within_trunc,
        "beyond_trunc": beyond_trunc,
        "max_abs": max_abs,
        "max_abs_step": max_abs_step,
        "max_rel_pct": max_rel,
        "max_rel_step": max_rel_step,
        "passed": passed,
        "ref_only": len(ref_q) - n,
        "test_only": len(test_q) - n,
    }
    return passed, stats


def print_stats(stats, label="", threshold=1e-6):
    if "error" in stats:
        print(f"  ERROR: {stats['error']}")
        return
    print(f"  Timesteps compared:   {stats['total']}")
    if stats["ref_only"] or stats["test_only"]:
        print(f"  Non-overlapping:      {stats['ref_only']} ref-only, {stats['test_only']} test-only")
    pct = lambda n: n / stats["total"] * 100 if stats["total"] else 0
    print(f"  Exact match (0 diff): {stats['exact']:6d}  ({pct(stats['exact']):6.2f}%)")
    print(f"  Within trunc noise:   {stats['within_trunc']:6d}  ({pct(stats['within_trunc']):6.2f}%)")
    print(f"  Beyond truncation:    {stats['beyond_trunc']:6d}  ({pct(stats['beyond_trunc']):6.2f}%)")
    print(f"  Max abs diff:         {stats['max_abs']:.2e} at step {stats['max_abs_step']}")
    print(f"  Max rel diff:         {stats['max_rel_pct']:.4f}% at step {stats['max_rel_step']}")
    status = "PASS" if stats["passed"] else "FAIL"
    print(f"  Result:               {status} (threshold: max abs diff < {threshold:.0e})")


def main():
    parser = argparse.ArgumentParser(description="Compare CFE discharge time series")
    parser.add_argument("reference", help="Reference discharge file")
    parser.add_argument("test_output", help="Test discharge file")
    parser.add_argument("--tol", type=float, default=1e-6, help="Absolute tolerance (default: 1e-6 m)")
    parser.add_argument("--trunc", type=float, default=5e-9, help="Truncation noise floor (default: 5e-9)")
    parser.add_argument("--ref-fmt", choices=["auto", "ngen", "standalone", "fred"], default="auto")
    parser.add_argument("--test-fmt", choices=["auto", "ngen", "standalone", "fred"], default="auto")
    args = parser.parse_args()

    for path, label in [(args.reference, "reference"), (args.test_output, "test")]:
        if not os.path.isfile(path):
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            return 2

    ref_fmt = args.ref_fmt if args.ref_fmt != "auto" else detect_format(args.reference)
    test_fmt = args.test_fmt if args.test_fmt != "auto" else detect_format(args.test_output)

    ref_q = READERS[ref_fmt](args.reference)
    test_q = READERS[test_fmt](args.test_output)

    label = os.path.basename(args.test_output)
    print(f"Comparing: {os.path.basename(args.reference)} ({ref_fmt}) vs {label} ({test_fmt})")

    passed, stats = compare(ref_q, test_q, args.tol, args.trunc)
    print_stats(stats, label, args.tol)

    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
