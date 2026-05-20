# samjr testmore status

Each subdirectory has a `script.R` derived from a corresponding SAM
`testmore` example, with `library(stockassessment)` swapped for
`library(samjr)`. The script runs the example, writes `res.out`, and the
test passes if `res.out` matches `res.EXP` byte-for-byte.

Run all scripts and tally results with:

```
./run-all.sh
```

or, in parallel:

```
make -j N testmore
```

from the repository root.

## Current status: 19 pass bit-identical to reference

`nscod`, `nscodFidx`, `nscodsw`, `nscodcovar`, `nscodXtraSd`, `nsher`,
`bfte2014`, `codIN3Ben`, `jit`, `mack`, `neaHaddockPredVar`,
`newtable`, `nscodbiopro`, `nscodswcwmofor`, `parallel`, `procres`,
`reduced`, `residplot`, `residuals`.

`newtable` compares samjr's tables against a snapshot of SAM
(`stockassessment 0.12.0`) tables saved as `newtable/tabSAM.rds`. The
test no longer needs SAM installed at run time; regenerate the snapshot
with `tools/build-tabSAM.R` or by running `extractTables(fitSAM)`
inside an interactive SAM session.

`nscodbiopro` and `nscodswcwmofor` use a samjr-generated `res.EXP`
rather than SAM's reference output. samjr's forecast uses the standard
year-y-1 survival convention `Z[y-1] = F[y-1] + M[y-1]` (consistent
with samjr's and SAM's own likelihoods), whereas SAM's `forecast.R`
mixes `Z = F[y-1] + nm[y]`. The discrepancy surfaces only when the
forecast horizon falls outside the natMor data range and uses
`ave.nm`.
