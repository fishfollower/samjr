## simulate-matidx.R
##
## Documents how IBTS_Q3_mat in survey.dat was generated. The series is
## derived from IBTS_Q3_gam by multiplying each (year, age) cell with the
## proportion mature MO[year, age], then perturbing with small log-normal
## noise. Running this script overwrites survey.dat, appending the
## simulated block after the two original surveys. Re-running with the
## same seed reproduces the values byte-equivalent (up to formatting).

library(samjr)

surveys <- read.ices("survey.dat")
mo      <- read.ices("mo.dat")

set.seed(123)
q3  <- surveys$IBTS_Q3_gam
yrs <- rownames(q3)
ags <- colnames(q3)
mox <- mo[yrs, ags]
matIdx <- q3 * mox * exp(matrix(rnorm(length(q3), sd = 0.1),
                                nrow = nrow(q3), ncol = ncol(q3)))

## Append a new survey block to survey.dat, preserving the two existing
## blocks byte-for-byte.
orig <- readLines("survey.dat")
tm   <- attr(q3, "time")
block <- c(
  "IBTS_Q3_mat",
  paste0(" ", min(as.integer(yrs)), " ", max(as.integer(yrs)), " "),
  paste0(" 1 1 ", tm[1], " ", tm[2], " "),
  paste0(" ", min(as.integer(ags)), " ", max(as.integer(ags)), " "),
  apply(matIdx, 1, function(r)
    paste(c("1", formatC(r, digits = 4, format = "f")), collapse = " "))
)
writeLines(c(orig, block), "survey.dat")
