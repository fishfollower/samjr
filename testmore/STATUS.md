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

## `ckmr`: regression-only, not a SAM comparison

`ckmr` exercises the close-kin mark-recapture feature: it fits nscod,
simulates a CKMR pair table from the fitted model with a fixed seed,
and refits with the parent-offspring term, the half-sibling term, and
both. Its `res.EXP` was generated from a verified samjr run - SAM has
no CKMR, so there is no external reference to compare against and this
test only guards against regression.

The `base` line reproduces `nscod`'s reference objective, `npar`
asserts that CKMR adds no estimated parameters (`psi` and the
abundance scale are fixed `conf` entries), and `ciratio` asserts that
the close-kin data sharpens terminal SSB.

The probabilities themselves were verified independently, against an
individual-based pedigree simulation rather than against SAM. The
half-sibling probability reproduces `ckmr/babyM.R`. The
parent-offspring probability agrees with `babyM.R` except when the
candidate parent is in the plus group, where `babyM.R` is biased low
(it back-calculates an exact age for a fish that could be older) and
can return a structural zero for pairs that genuinely occur; samjr
averages fecundity over the possible ages instead. See
`samjr/R/ckmr.R`.

## `babyM`: the close-kin runs from `ckmr/babyM.R`

`babyM` reproduces the three runs at the end of `ckmr/babyM.R` - official data
plus POP, plus HSP, plus both - on the mackerel assessment
(`model_XplConf.08`) that script was built around. The data was extracted from
that fit into ordinary ICES `.dat` files, so the case is self-contained and
needs neither the 6.7 MB `RData` nor the `stockassessment` package.

The `base` line is a genuine external check: it reproduces SAM's own objective
for that assessment, 580.2740544, to within 2e-8 with the same 40 parameters
and the same 26 imputed observations. The close-kin lines are
regression-only, as for `ckmr`.

The case writes `babyM/babyM.png`: three panels of SSB, the official
assessment band in red with the close-kin fit overlaid in blue and the
estimated `psi` in the corner, following the three panels `babyM.R` draws at
its end. The legend keys on the confidence bands rather than the lines,
because what the close-kin data buys is width rather than level. Like
`res.out` it is a generated artefact, not committed.

### Three documented differences from the prototype

* `babyM.R` also estimates a global natural-mortality multiplier in its two
  HSP runs. samjr has no such parameter, so those runs differ by that one
  degree of freedom. `psi` is estimated in all three, as in `babyM.R`
  (`conf$ckmrEstimatePsi`).

* The runs are seeded from the baseline assessment rather than from `defpar`
  defaults. Without the mortality multiplier the HSP likelihood is bimodal in
  `psi`, and a cold start converges to a psi -> infinity mode with a worse
  objective (785.3 against 741.4). This is a genuine property of the
  likelihood, not a solver failure - the bad mode has a max absolute gradient
  of 4e-8.

* The genotyped sample is uniform over ages, where `babyM.R` took ages 2-4
  only (`selection=c(1,1,1,0,...)`). That window is a poor design for
  parent-offspring pairs and it never exercises the plus-group code path; see
  below.

### Why the sampling design was widened

Under `babyM.R`'s three-age window the POP-only panel looks as though the
close-kin data does nothing. It is not a bug - POP alone narrows the mean SSB
interval to 0.984 of the baseline, against 0.945 for HSP alone - but two
things hold it back, both properties of that window:

* `P_POP` is a fecundity ratio over an abundance, `2 fec(a_p, b) / TRO(b)`.
  Sampling only ages 2-4 puts every candidate parent on a nearly flat stretch
  of the fecundity schedule, so `psi` and abundance trade off almost
  perfectly. Holding `psi` fixed takes POP from 0.984 to 0.870.
* Only 30 of the 66 pairs contribute to POP at all, with parent ages at the
  birth year spanning 2-4 and birth years 2018-2021.

Sampling uniformly removes both effects: 746 contributing pairs, parent ages
0-12, birth years 2009-2021, and POP alone reaches 0.819 whether `psi` is free
(0.819) or fixed (0.822). POP and HSP together still beat either alone or the
product of the two, because the half-sibling term identifies `psi` and so
frees the POP term to inform abundance.

Uniform sampling also puts candidate parents in the plus group, which the
narrow window never does. Of 946 pairs, 158 need the plus-group POP treatment
and 105 of the 383 observed POPs have a plus-group parent. Of those, 19 pairs
- carrying 5 actually observed POPs - would have had probability exactly zero
under `babyM.R`'s back-calculated age, each contributing about 230 to the
objective through its `1e-100` guard.

### Choice of seed

The estimates are unbiased. Over 8 seeds the mean SSB shift relative to the
baseline is +0.9% for POP, +0.4% for HSP and +1.0% for both, none
distinguishable from zero, and `psi` averages 1.49 to 1.51 against a simulated
1.5. Single realisations move around a fair amount, though: over 20 seeds the
peak-year SSB shift spans -7.4% to +10.5% for POP and -4.1% to +4.4% for HSP,
and `psi` spans 1.11 to 1.96.

The seed is therefore chosen from a 20-seed scan as the draw closest to the
median on all six of those statistics - the most typical draw, not the most
favourable. An earlier seed happened to sit 4.1 sd below the mean on the
HSP-only SSB shift (-9.5%, outside the range of all 20 others) and produced a
figure suggesting the half-sibling data moves the assessment about twice as
much as it typically does. That was not a convergence artefact: refitting that
data set from eight jittered starts reached the identical optimum every time.

## `babyMibm`: the same runs, checked against a real pedigree

`babyMibm` is `babyM` with `simulateCKMR` replaced by an individual-based
pedigree simulation in `babyMibm/ibmsim.R`, local to that folder. The IBM
builds explicit individuals with recorded mothers and fathers and counts the
kin pairs that actually exist among a sample of them. It never evaluates a
close-kin probability, so comparing its kin counts against `ckmrProb` is an
independent check of the formulas - something `simulateCKMR` cannot give,
since it draws from those same formulas and is self-consistent by
construction.

Two things had to be got right for the comparison to mean anything:

* **The population must satisfy the survival recursion.** SAM's fitted `N`
  carries process noise, so realised survival ratios
  `N[y+1,a+1]/(N[y,a]exp(-Z))` run from 0.8 to 5.6 and compound to four orders
  of magnitude along a cohort track. The close-kin probabilities mix `N` and
  `Z`, so an IBM built on the fitted `N` would be testing that inconsistency
  rather than the formulas. `ibmConsistentN` rebuilds `N` from the fitted
  recruitment by deterministic survival; SSB moves by at most 9%.
* **Sampling must be lethal.** Drawing each sample year independently lets the
  same fish be caught twice and appear at two ages, which contradicts the
  assumption the probabilities rest on (an adult cannot parent anything born
  after it was sampled) and undercounts kin. 4.5% of sample rows were repeats,
  and removing them moved the POP ratio from 0.937 to 0.997.

The true stock is about 8e11 individuals, so the population is scaled down by
1e5 and `conf$ckmrScale` set to 0.01 to match; model numbers-at-age times
`ckmrScale` are then exactly the simulated individuals. Full siblings, which
the factor-4 half-sibling expression assumes away, do not appear even at that
scale.

### What the check finds

Where ages are known exactly the formulas reproduce the pedigree. Over ten
independent pedigrees the ratio of observed to predicted kin is 0.997
(p = 0.72) for parent-offspring pairs and 0.999 (p = 0.90) for half-siblings
on rows with no plus-group animal. The plus-group POP treatment also holds up,
at 1.009.

Half-siblings involving a plus-group animal are a different matter: the ratio
there is about 0.75. This is a real limitation, not a defect in the
implementation. 58% of sampled age-12 fish are in fact older, because the plus
group lumps every age from 12 upward, so their recorded birth year is too late
and the model computes the wrong cohort gap. The POP expression averages over
the ages a plus-group *parent* could have had, which is why it survives; there
is no equivalent averaging over a plus-group animal's own birth year in the
half-sibling term, and adding one would require carrying a distribution over
its true age through both `b1` and `b2`. Until then, close-kin sampling
designs should treat plus-group animals with caution, or exclude them.

The IBM also produces a small number of parent-offspring pairs within a single
(year, age) cell - impossible for known ages, but possible in the plus group,
where two fish recorded at the same age can differ in real age by enough for
one to have parented the other. The pair table has no row for a cell paired
with itself, so those are counted and reported separately rather than used.

`res.EXP` records the observed-to-predicted ratios alongside the three run
objectives, and asserts `ok.pop` and `ok.hsp`: the formulas must reproduce the
pedigree to within 10% on rows where ages are known.
