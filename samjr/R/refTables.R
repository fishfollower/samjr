## Print / summary methods for sam_referencepoints.

##' @method print sam_referencepoints
##' @export
print.sam_referencepoints <- function(x, ...){
  cat("samref reference points",
      if(isTRUE(x$stochastic)) "(stochastic)" else "(deterministic)", "\n")
  cat("Catch type:", x$catchType, " aveYears:",
      paste(range(x$aveYears), collapse = "-"),
      " selYears:", paste(range(x$selYears), collapse = "-"), "\n\n")
  if(length(x$tables) == 0L || nrow(x$tables$F) == 0L){
    cat("(no reference points)\n"); return(invisible(x))
  }
  print(round(x$tables$F, 4))
  invisible(x)
}

##' @method summary sam_referencepoints
##' @export
summary.sam_referencepoints <- function(object, ...){
  blocks <- lapply(names(object$tables), function(nm){
    m <- object$tables[[nm]]
    rownames(m) <- paste(rownames(m), nm, sep = ".")
    m
  })
  do.call(rbind, blocks)
}
