# samjr testmore status

Each subdirectory has a `script.R`, most of them derived from a
corresponding SAM `testmore` example with `library(stockassessment)`
swapped for `library(samjr)`. The script runs the example, writes
`res.out`, and the test passes if `res.out` matches `res.EXP`
byte-for-byte.

Run all scripts and tally results with:

```
./run-all.sh
```

or, in parallel:

```
make -j N testmore
```

from the repository root.

## Current status: all 33 pass

As of 2026-11-21 (`make testmore`).

### Fits, tables, residuals and forecasts

`nscod`, `nscodFidx`, `nscodsw`, `nscodcovar`, `nscodXtraSd`,
`nscodMatIdx`, `nscodplusgroup`, `nsher`, `bfte2014`, `codIN3Ben`,
`jit`, `mack`, `neaHaddockPredVar`, `newtable`, `nscodbiopro`,
`nscodswcwmofor`, `parallel`, `procres`, `reduced`, `residplot`,
`residuals`, `whbrecpred`.

`nscodMatIdx` adds a mature-fish index (fleet type 4). `nscodplusgroup`
fits with the plus group switched off for one fleet
(`conf$maxAgePlusGroup <- c(1, 1, 0)`).

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

`whbrecpred` (blue whiting) is samjr-only: it fits recruitment codes 0,
62 and 102 (wsc scaled by 1e6), peels 10 years with `retro`, predicts
recruitment one and two years ahead from each peel with `forecast(...,
useRproc = TRUE)`, and writes the RMSE of log(predicted / full run) per
code and horizon. Runs in about 65s. `plotscript.R` (not part of the
test) overrides `srCodes`, `nPeel` and `resFile` before sourcing
`script.R` to run codes 0, 1, 2, 62 and 102 with 15 peels, writes the
RMSEs to `recpred.out` and plots both horizons to `recpred.pdf` (about
150s).

### Reference points and harvest control rules

`cis`, `crash_ext`, `msy`, `perrecruit`, `xspr_xb0` and `yprparity`
compare `referencepoints` / `deterministicReferencepoints` estimates
with SAM 0.12.0 output; `cis` checks the confidence intervals within a
tolerance. `hcr` checks the trapezoidal rule against exact algebraic
answers, and `hcr_plots` checks that the HCR plot methods run without
error.

`hcrsam` compares `hcrFun` with SAM's rule function (`C_hcrR`) on a grid
of SSB values for four rule settings, including `Forigin > 0`. SAM's
values come from the snapshot `hcrSAM.rds`; rebuild it with
`make testmore/hcrsam/hcrSAM.rds` (needs stockassessment). Only the rule
is compared, because SAM 0.12.0's `hcr()` projection stops with
"Forecast type not implemented".

`hcrcheck` checks the `hcr()` projection with recruitment noise switched
off: at constant F it settles at the per-recruit SPR and YPR (nscod
random-walk fit) and at the Ricker and Beverton-Holt equilibria (NEA
saithe, data from `customsr`), and each year's Fbar equals
`hcrFun(SSB)` under the rule.

### Custom recruitment

`customsr` checks that a user-supplied `logSRfun` / `logSRinit`
reproduces the built-in Ricker and Beverton-Holt fits (NEA saithe,
1960-2020): same nll and same SR parameters at convergence.
