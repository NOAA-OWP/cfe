#!/bin/bash
# compare_outputs.sh — numerically compare two CFE output files
# Usage: compare_outputs.sh <reference> <test_output> [tolerance] [label]
# Returns: 0 on PASS, 1 on FAIL, 2 on input error
set -e
REF="$1"; TST="$2"; TOL="${3:-1e-7}"; LABEL="${4:-$(basename "$TST")}"
[ -z "$REF" ] || [ -z "$TST" ] && { echo "Usage: compare_outputs.sh <ref> <test> [tol] [label]" >&2; exit 2; }
[ -f "$REF" ] || { echo "ERROR: reference not found: $REF" >&2; exit 2; }
[ -f "$TST" ] || { echo "ERROR: test output not found: $TST" >&2; exit 2; }
REF_CLEAN=$(mktemp); TST_CLEAN=$(mktemp); trap 'rm -f "$REF_CLEAN" "$TST_CLEAN"' EXIT
grep -v '^[[:space:]]*#' "$REF" | grep -v '^[[:space:]]*$' > "$REF_CLEAN"
grep -v '^[[:space:]]*#' "$TST" | grep -v '^[[:space:]]*$' > "$TST_CLEAN"
awk -v tol="$TOL" -v label="$LABEL" '
BEGIN { max_err=0; n_vals=0; n_fail=0 }
{
    if ((getline tline < tst_file) <= 0) { print "ERROR: test file short at line " NR; exit 3 }
    n=split($0,rfields); m=split(tline,tfields)
    if (n!=m) { printf "FAIL [%s]: field count mismatch line %d: ref=%d test=%d\n",label,NR,n,m; exit 1 }
    for (i=1;i<=n;i++) {
        rv=rfields[i]+0; tv=tfields[i]+0; rv_abs=rv<0?-rv:rv
        if (rv_abs<1e-15) err=tv<0?-tv:tv
        else { diff=rv-tv; if(diff<0)diff=-diff; err=diff/rv_abs }
        n_vals++; if(err>max_err)max_err=err
        if (err>tol) { n_fail++
            if(n_fail<=5) printf "  line %d, field %d: ref=%.8e got=%.8e (rel_err=%.2e)\n",NR,i,rv,tv,err
            if(n_fail==6) print "  (further mismatches suppressed)"
        }
    }
}
END {
    if (n_fail>0) { printf "FAIL [%s]: %d/%d exceed tol=%.1e, max_err=%.2e\n",label,n_fail,n_vals,tol,max_err; exit 1 }
    else { printf "PASS [%s]: %d values, max_rel_err=%.2e (tol=%.1e)\n",label,n_vals,max_err,tol }
}' tst_file="$TST_CLEAN" "$REF_CLEAN"
