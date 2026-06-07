# samjr

A minimal, pure-R reimplementation of the
[State-space Assessment Model (SAM)](https://github.com/fishfollower/SAM)
built on top of [RTMB](https://kaskr.r-universe.dev/RTMB).
samjr trades C++ for readable R while keeping bit-identical fits to SAM
on every test case shipped with the package.

> **Status.** A small subset of SAM, with the canonical North Sea cod
> workflow as the reference example. All 19 testmore scripts pass
> bit-identical to a SAM 0.12.0 reference.

## Why samjr?

* **Readable.** The full state-space likelihood lives in
  [`samjr/R/likelihood.R`](samjr/R/likelihood.R) (~350 lines of base R).
* **Hackable.** Want to tweak the recruitment model, the F correlation,
  or add a new biology process? Edit one R file and reinstall. No
  recompilation, no TMB template gymnastics.
* **Faithful.** The default workflow (`read.ices` → `setup.sam.data` →
  `defcon` → `defpar` → `sam.fit`) and the standard tables (`ssbtable`,
  `fbartable`, `rectable`, `catchtable`, `tsbtable`, `ntable`, `faytable`,
  `caytable`, `qtable`, `partable`, `modeltable`) all match SAM to within
  `tolerance = 1e-4`.
* **RTMB-native.** Random effects, Laplace approximation, OSA residuals,
  joint precision, automatic differentiation - all from RTMB.

## Installation

```r
# from a local clone:
install.packages(".", repos = NULL, type = "source")
# or:
remotes::install_github("fishfollower/samjr")
```

samjr depends on `RTMB` and uses base R for everything else (`Matrix`,
`methods`, `stats`, `utils`, `graphics`, `grDevices`, `parallel`).

## Quick start

The package ships with the North Sea cod (`nscod`) example baked in:

```r
library(samjr)
data(nscodData)
data(nscodConf)
par <- defpar(nscodData, nscodConf)
fit <- sam.fit(nscodData, nscodConf, par)

fit
#> samjr model: log likelihood is -145.5167  Convergence OK

opar <- par(mfrow = c(2, 2))
ssbplot(fit); fbarplot(fit); recplot(fit); catchplot(fit)
par(opar)

ssbtable(fit)        # Estimate / Low / High by year
fbartable(fit)
catchtable(fit)
partable(fit)        # parameter estimates with std. errors
modeltable(fit)      # one-line model summary (nll, npar, AIC)

res <- residuals(fit)  # OSA residuals via RTMB::oneStepPredict
plot(res)
```

To start from your own ICES `.dat` files instead:

```r
cn  <- read.ices("cn.dat");   cw  <- read.ices("cw.dat")
sw  <- read.ices("sw.dat");   nm  <- read.ices("nm.dat")
mo  <- read.ices("mo.dat");   pf  <- read.ices("pf.dat")
pm  <- read.ices("pm.dat")
surveys <- read.ices("survey.dat")

dat  <- setup.sam.data(surveys = surveys, residual.fleet = cn,
                       prop.mature = mo, stock.mean.weight = sw,
                       catch.mean.weight = cw, prop.f = pf, prop.m = pm,
                       natural.mortality = nm)
conf <- defcon(dat)
par  <- defpar(dat, conf)
fit  <- sam.fit(dat, conf, par)
```

A short stochastic forecast:

```r
set.seed(123)
fc <- forecast(fit, fscale = c(1, 1, 1, 1))
fc
```

## The model

samjr is a state-space age-structured assessment model. Latent states
are log numbers-at-age $\log N_{a,y}$ and log fishing mortality
$\log F_{a,y}$; observations are catches-at-age and survey indices.

### Population dynamics

Numbers at age evolve as

```math
\log N_{a,y} =
\begin{cases}
  R(\mathrm{SSB}_{y-a}) + \varepsilon^{R}_y & a = a_{\min} \\
  \log N_{a-1,\, y-1} - F_{a-1,\, y-1} - M_{a-1,\, y-1} + \varepsilon^{S}_{a,y}
    & a_{\min} < a < A \\
  \log(N_{A-1,\,y-1}\,e^{-Z_{A-1,y-1}} + N_{A,\,y-1}\,e^{-Z_{A,y-1}}) + \varepsilon^{S}_{A,y}
    & a = A \quad (\text{plus group})
\end{cases}
```

with $Z_{a,y} = F_{a,y} + M_{a,y}$, and $a_{\min}$ the first modelled
age (set by `conf$minAge`, not necessarily 1). The recruitment function
$R(\cdot)$ is one of:

* `srmode = 0`: random walk, $R(\cdot) = \log N_{a_{\min},\,y-1}$;
* `srmode = 1`: Ricker, $R(S) = \alpha + \log S - e^{\beta} S$;
* `srmode = 2`: Beverton-Holt, $R(S) = \alpha + \log S - \log(1 + e^{\beta} S)$.

Process noise:

```math
\varepsilon^{R}_y \sim \mathcal{N}(0, \sigma_R^2),
\qquad
\varepsilon^{S}_{a,y} \sim \mathcal{N}(0, \sigma_S^2).
```

Fishing mortality is a multivariate random walk in log space,

```math
\log F_y \;\sim\; \mathcal{N}(\log F_{y-1},\, \Sigma_F),
\qquad
[\Sigma_F]_{ij} = \sigma_{F,i}\,\sigma_{F,j}\,\rho_{ij},
```

with three options for the correlation structure
$\rho_{ij}$ (`fcormode`):

* `0` independent: $\rho_{ij} = \mathbf{1}\{i = j\}$;
* `1` compound symmetry: $\rho_{ij} = \rho$;
* `2` AR(1): $\rho_{ij} = \rho^{|i-j|}$.

### Observations

Observation type is set per fleet via `fleetType` and modelled on the
log scale. Stacking the log observations for fleet $f$ in year $y$ into
a vector $\log \widetilde O_{f,y} \in \mathbb{R}^{n_{f,y}}$ (one entry
per sampled age), the observation model is

```math
\log \widetilde O_{f,y} \;\sim\; \mathcal{N}(\log O_{f,y},\, \Sigma_{f,y}),
```

where $\log O_{f,y} \in \mathbb{R}^{n_{f,y}}$ is the vector of predicted
log observations and $\Sigma_{f,y}$ is $n_{f,y} \times n_{f,y}$. The
scalar prediction at age $a$ is

```math
[\log O_{f,y}]_a =
\begin{cases}
  \log N_{a,y} + \log F_{a,y} - \log Z_{a,y} + \log(1 - e^{-Z_{a,y}})
    & \text{catch (type 0)} \\
  \log Q_{a,f} + \log N_{a,y} - \tau_f Z_{a,y}
    & \text{survey index (type 2)} \\
  \log Q_f + \log \mathrm{SSB}_y
    & \text{biomass index (type 3)} \\
  \log Q_{a,f} + \log \mathrm{MO}_{a,y} + \log N_{a,y} - \tau_f Z_{a,y}
    & \text{mature-fish index (type 4)}
\end{cases}
```

where $\tau_f$ is the survey sample time and $Q$ is the catchability.
Biomass indices (type 3) carry no age dimension, so $\log O_{f,y}$ and
$\Sigma_{f,y}$ collapse to scalars (i.e. $n_{f,y} = 1$).

The covariance factorises as $\Sigma_{f,y} = D_{f,y} R_{f,y} D_{f,y}$,
with correlation $R_{f,y}$ (independent, IGAR-distance, or unstructured
per fleet via `conf$obsCorStruct`) and a per-observation standard
deviation on the diagonal. The SD defaults to
$\exp([\texttt{logSdLogObs}]_{\texttt{keyVarObs}})$ but, as in SAM, can be
overridden per observation by a supplied weight (`conf$fixVarToWeight`,
from a `weight`/`cov`/`cov-weight` attribute), an extra-SD factor
(`conf$keyXtraSd`), or a prediction-variance link that grows the SD with
the prediction (`conf$predVarObsLink`). These modifiers only set the
diagonal, so they compose with any correlation structure.

### Spawning stock biomass and catch in weight

```math
\mathrm{SSB}_y = \sum_a N_{a,y}\, \mathrm{SW}_{a,y}\, \mathrm{MO}_{a,y}
\,e^{-\mathrm{PF}_{a,y} F_{a,y} - \mathrm{PM}_{a,y} M_{a,y}}
```

```math
\mathrm{Catch}_y = \sum_a N_{a,y}\,\mathrm{CW}_{a,y}\,
\frac{F_{a,y}}{Z_{a,y}}(1 - e^{-Z_{a,y}})
\qquad
\overline{F}_y = \frac{1}{|A_{\bar F}|} \sum_{a \in A_{\bar F}} F_{a,y}
```

### Optional GMRF biology processes

Stock weight $\mathrm{SW}$, catch weight $\mathrm{CW}$, maturity
$\mathrm{MO}$, and natural mortality $M$ can each be modelled as
year-by-age Gaussian Markov random fields rather than treated as
data. The latent matrix $X$ then has process precision

```math
Q = I - \phi_{\text{cohort}} W_c - \phi_{\text{within-year}} W_d
\quad (- \phi_{\text{between-year}} W_p),
```

with sparse band matrices $W_c, W_d, W_p$ encoding cohort, within-year,
and between-year neighbours. Observations on positive-valued biology
(SW, CW, M) are log-normal; maturity uses a beta likelihood with the
SAM "squash" trick so $0,1$ proportions stay strictly inside $(0,1)$.
Toggle each process with `conf$stockWeightModel`,
`conf$catchWeightModel`, `conf$matureModel`,
`conf$mortalityModel` (`0` = data, `1` or `2` = GMRF).

### Estimation

Random effects are
$\log N$, $\log F$, plus the GMRF latent matrices and any imputed
missing observations. Fixed effects are the sd / correlation
parameters, $Q$, recruitment hyperparameters, and the GMRF means and
$\phi$'s. Estimation is the standard Laplace approximation provided by
RTMB:

```r
obj <- RTMB::MakeADFun(makeLikelihood(dat), parameters,
                       random = c("logN", "logF", ...))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- RTMB::sdreport(obj)
```

`sam.fit` wraps that with sensible defaults (3 Newton steps,
parameter-name-aware bounds, and an `sdreport` with joint precision).

## Tables, plots, and forecasting

| function | what it returns |
| --- | --- |
| `ssbtable`, `fbartable`, `rectable`, `catchtable`, `tsbtable` | year x (Estimate, Low, High) |
| `ntable`, `faytable`, `caytable` | year x age, point estimates |
| `qtable` | catchabilities by fleet x age |
| `partable`, `modeltable` | parameter / model summaries |
| `modelDescription` | prose description of the fitted model |
| `ssbplot`, `fbarplot`, `recplot`, `catchplot`, `tsbplot` | summary plots |
| `selplot`, `srplot`, `fitplot`, `dataplot`, `parplot`, `sdplot`, `corplot` | diagnostic plots |
| `residuals` (and `plot`/`print` methods) | OSA residuals via `RTMB::oneStepPredict` |
| `procres` | joint-sample process residuals for $\log N$ and $\log F$ |
| `residplot` | p-value heatmaps (bias, variance, age/time correlation, mean-variance, normality) |
| `forecast` | short-term stochastic forecast |
| `ypr`, `yprtable`, `yprplot` | yield-per-recruit analysis |
| `retro`, `runwithout`, `leaveout`, `mohn` | retrospective and leave-one-out tools |
| `jit` | jitter starting values, refit |
| `simulate`, `simstudy` | simulate from the fitted model / parametric bootstrap |
| `c.sam` | combine multiple `sam` fits into a `samset` |
| `coef`, `logLik`, `nobs`, `AIC`, `summary` | standard S3 methods on a `sam` fit |
| `read.ices`, `read.surveys`, `read.data.files` | read ICES / SAM input files |
| `write.ices`, `write.surveys`, `write.data.files` | inverse writers |
| `loadConf`, `saveConf` | round-trip a `conf` list to / from disk |
| `getFleet`, `reduce` | extract a fleet / subset a dataset |

`forecast.sam` projects forward year by year using the standard
$\mathrm{SAM}$ survival recursion ($Z = F_{y-1} + M_{y-1}$). The
fitted joint precision is used to draw initial state perturbations,
recruitment is resampled from `rec.years`, and biology can be either
averaged over `ave.years` or sampled from the GMRF process models.

## Repository layout

```
samjr/
  samjr/             # the R package
    R/               # source (see R/likelihood.R for the model)
    data/            # nscodData, nscodConf, nscodParameters
    man/             # roxygen-generated .Rd files
    tools/           # build-nscod-data.R
  testmore/          # 19 end-to-end test cases (script.R + res.EXP)
  Makefile           # see `make help`
```

## Development

```sh
make doc        # regenerate man/*.Rd and NAMESPACE from roxygen
make install    # doc + R CMD INSTALL
make build      # doc + R CMD build (produces samjr_*.tar.gz)
make check      # build + R CMD check --as-cran
make data       # rebuild data/nscod*.rda from testmore/nscod
make testmore   # run every testmore/<dir>/script.R; pass = bit-identical
```

`make -j N testmore` runs the testmore suite in parallel. Each test is
self-contained (no shared state), so parallel speedup is essentially
linear in cores.

## Acknowledgements

samjr exists because of the work that went into
[SAM](https://github.com/fishfollower/SAM) and
[RTMB](https://kaskr.r-universe.dev/RTMB) - this package is a small
pure-R restatement of the SAM model on top of RTMB and owes everything
to those projects.

## License

GPL-2.
