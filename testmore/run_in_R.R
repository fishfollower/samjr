
#remotes::install_github("fishfollower/samjr/samjr@dsem")

library(testthat)
library(samjr)

setwd( R'(C:\Users\jtuth\Documents\GitHub\fishfollower\samjr)' )

# Read a file as lines with any trailing CR removed (like diff --strip-trailing-cr)
read_lines_nocr <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
# Test if runs OK
runs_ok <- function(expr) {
  tryCatch({ expr; TRUE },
           error = function(e) structure(FALSE, msg = conditionMessage(e)))
}

# tests/testthat/test-testmore.R
dirs <- list.dirs("testmore", recursive = FALSE)
#dirs = dirs[-(1:17)]

# FAILS
#perrecruit
#reduced
tests = data.frame(
  dirs = dirs,
  runs = rep(NA, length(dirs)),
  matches = rep(NA, length(dirs))
)

for (i in seq_along(dirs)) {
  message("running ", dirs[i])
  tests$runs[i] = runs_ok(callr::r(
    function(dir){
      setwd(dir);
      source("script.R")
    },
    args = list(dirs[i]),
    stdout = file.path(dirs[i], "stdout.log"),
    stderr = file.path(dirs[i], "stdout.log")
  ))

  if( isTRUE(tests$runs[i]) ){
    res_out = read_lines_nocr(file.path(dirs[i], "res.out"))
    res_exp = read_lines_nocr(file.path(dirs[i], "res.EXP"))
    tests$matches[i] = identical( res_out, res_exp )
  }
}


