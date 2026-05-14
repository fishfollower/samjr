#!/bin/bash
# Run each OK testmore script under samjr and stockassessment, time each.
# Usage: bash bench_ok.sh
set -u
TESTS=(bfte2014 jacobian jit mack neaHaddockPredVar
       nscod nscodFidx nscodcovar nscodsw nscodXtraSd nsher reduced residuals)

ROOT="$(cd "$(dirname "$0")" && pwd)"
RESULT="$ROOT/bench_ok_result.txt"
: > "$RESULT"
printf "%-22s %12s %12s %8s\n" test samjr_s stockassessment_s ratio | tee -a "$RESULT"

run_one() {
  local dir="$1" pkg="$2"
  ( cd "$ROOT/$dir" && rm -f res.out
    /usr/bin/time -f "%e" -o /tmp/__bench_t.txt \
      Rscript --vanilla -e "
        suppressMessages(library($pkg))
        # The scripts call library(samjr) themselves; pre-load $pkg so they share NS
        src <- readLines('script.R')
        src <- sub('library\\\\(samjr\\\\)', 'library($pkg)', src)
        eval(parse(text=src))
      " >/dev/null 2>/dev/null
    cat /tmp/__bench_t.txt
  )
}

for t in "${TESTS[@]}"; do
  bb=$(run_one "$t" samjr)
  sa=$(run_one "$t" stockassessment)
  ratio=$(awk -v a="$bb" -v b="$sa" 'BEGIN{ if (b>0) printf "%.2fx", a/b; else print "n/a" }')
  printf "%-22s %12s %12s %8s\n" "$t" "$bb" "$sa" "$ratio" | tee -a "$RESULT"
done
