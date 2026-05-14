#!/bin/bash
# Run every samjr testmore example and tally results.
# Usage: ./run-all.sh
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT" || exit 1
PASS=(); FAIL_DIFF=(); FAIL_ERR=()
declare -A TIMES

elapsed() {
  awk -v s="$1" -v e="$2" 'BEGIN{printf "%.1f", e-s}'
}

for dir in */; do
  dir="${dir%/}"
  [ -f "$dir/script.R" ] || continue
  t0=$(date +%s.%N)
  ( cd "$dir" && rm -f res.out errlog.txt
    timeout 120 R --slave --vanilla -e 'source("script.R")' \
      > /dev/null 2> errlog.txt
  )
  t1=$(date +%s.%N)
  TIMES[$dir]=$(elapsed "$t0" "$t1")
  if [ -f "$dir/res.out" ] && [ -f "$dir/res.EXP" ]; then
    if diff --strip-trailing-cr -q "$dir/res.out" "$dir/res.EXP" > /dev/null; then
      PASS+=("$dir")
    else
      FAIL_DIFF+=("$dir")
    fi
  else
    FAIL_ERR+=("$dir")
  fi
done

print_with_time() {
  for d in "$@"; do
    printf '  %-25s %6ss\n' "$d" "${TIMES[$d]:-?}"
  done
}

echo "===== PASS - bit-identical to reference (${#PASS[@]}) ====="
print_with_time "${PASS[@]}"
echo
echo "===== RAN, RESULT DIFFERS (${#FAIL_DIFF[@]}) ====="
print_with_time "${FAIL_DIFF[@]}"
echo
echo "===== ERRORED (${#FAIL_ERR[@]}) ====="
for d in "${FAIL_ERR[@]}"; do
  err=$(grep -m1 -E '^Error' "$d/errlog.txt" 2>/dev/null | head -1 | cut -c1-100)
  printf '  %-25s %6ss  %s\n' "$d" "${TIMES[$d]:-?}" "${err:-no error captured}"
done
echo

total_time=0
for d in "${!TIMES[@]}"; do
  total_time=$(awk -v a="$total_time" -v b="${TIMES[$d]}" 'BEGIN{printf "%.1f", a+b}')
done

echo "===== TOTAL ====="
echo "  scripts run     : $((${#PASS[@]} + ${#FAIL_DIFF[@]} + ${#FAIL_ERR[@]}))"
echo "  passed          : ${#PASS[@]}"
echo "  ran, diff fail  : ${#FAIL_DIFF[@]}"
echo "  errored         : ${#FAIL_ERR[@]}"
echo "  wall-clock total: ${total_time}s"
