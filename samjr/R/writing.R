##' Write a matrix to an ICES/CEFAS data file
##'
##' Takes a year-by-age matrix (with year row names and age column names)
##' and writes it in the ICES/CEFAS \code{.dat} format. Custom attributes
##' on \code{x} (anything other than \code{dim}/\code{dimnames}) are
##' appended as \code{#'@name} / \code{#' value} comment lines so they
##' round-trip through \code{\link{read.ices}}.
##'
##' If \code{x} is a length-1 list of matrices, the single matrix is
##' written. If \code{x} is a longer list (or a 3D array), each matrix is
##' written to a separate file using a \code{_NNNNN} suffix on the file
##' stem - this is the SAM convention; samjr itself only ships one
##' residual catch fleet but the writer handles either form so a
##' setup-and-write round-trip stays compatible.
##'
##' @param x a matrix, a 3D array, or a list of matrices.
##' @param fileout output file path.
##' @param writeToOne if \code{TRUE} and all matrices in a list/array are
##' equal, only one file is written.
##' @param ... arguments forwarded to \code{\link{write}}.
##' @return invisibly \code{NULL}.
##' @export
write.ices <- function(x, fileout, writeToOne = TRUE, ...){
  writeOne <- function(x, fileout, ...){
    top <- paste0(fileout, " auto written\n1 2\n",
                  paste0(range(as.integer(rownames(x))), collapse = "  "), "\n",
                  paste0(range(as.integer(colnames(x))), collapse = "  "), "\n1")
    write(top, fileout, ...)
    write(t(x), fileout, ncolumns = ncol(x), append = TRUE, sep = "  \t", ...)
    extra <- !(names(attributes(x)) %in% c("dim", "dimnames"))
    if(any(extra)){
      a <- attributes(x)[extra]
      for(i in seq_along(a)){
        write(paste0(sprintf("#'@%s\n", names(a)[i]),
                     paste(paste("##'", deparse(a[[i]])), collapse = "\n")),
              fileout, append = TRUE, ...)
      }
    }
  }
  if(is.matrix(x)){
    writeOne(x, fileout, ...)
    return(invisible(NULL))
  }
  if(is.array(x) && length(dim(x)) == 3){
    x <- lapply(split(x, rep(seq_len(dim(x)[3]), each = prod(dim(x)[1:2]))),
                matrix, nrow = dim(x)[1], ncol = dim(x)[2],
                dimnames = dimnames(x)[1:2])
  }
  if(!is.list(x))
    stop("x must be a matrix, an array or a list of matrices")
  if(writeToOne && all(vapply(x, function(y) isTRUE(all.equal(y, x[[1]])),
                              logical(1))))
    x <- x[1]
  if(length(x) == 1){
    writeOne(x[[1]], fileout)
  }else{
    f2 <- sub("\\.dat", "_%05d.dat", fileout)
    for(i in seq_along(x)) writeOne(x[[i]], sprintf(f2, i))
  }
  invisible(NULL)
}

##' Extract the observed (year x age) matrix for a single fleet
##' @param data a \code{sam_data} list (or a fitted \code{sam} object;
##' the \code{$data} slot is used in that case).
##' @param fleet integer fleet index.
##' @return a matrix with year row names and age column names; cells
##' without an observation are zero.
##' @export
getFleet <- function(data, fleet){
  if(inherits(data, "sam")) data <- data$data
  fidx <- data$aux[, "fleet"] == fleet
  aux  <- data$aux[fidx, , drop = FALSE]
  logo <- data$logobs[fidx]
  goget <- function(y, a){
    ret <- exp(logo[aux[, "year"] == y & aux[, "age"] == a])
    if(length(ret) == 0) 0 else ret
  }
  yr <- min(aux[, "year"]):max(aux[, "year"])
  ar <- min(aux[, "age"]):max(aux[, "age"])
  tmp <- outer(yr, ar, Vectorize(goget))
  dimnames(tmp) <- list(yr, ar)
  tmp
}

##' Write the survey fleets of a samjr data object to an ICES survey file
##' @param data a \code{sam_data} list or a fitted \code{sam} object.
##' @param fileout output file path.
##' @param ... arguments forwarded to \code{\link{write}}.
##' @return invisibly \code{NULL}.
##' @export
write.surveys <- function(data, fileout, ...){
  if(inherits(data, "sam")) data <- data$data
  sidx <- which(data$fleetTypes %in% c(2, 3, 4))
  top <- paste0(fileout, " auto written\n", 100 + length(sidx))
  write(top, fileout, ...)
  fleetNames <- attr(data, "fleetNames")
  for(s in sidx){
    write(fleetNames[s], fileout, append = TRUE, ...)
    S  <- getFleet(data, s)
    yr <- range(as.integer(rownames(S)))
    ar <- range(as.integer(colnames(S)))
    write(paste0(yr[1], " ", yr[2]), fileout, append = TRUE, ...)
    st <- data$sampleTimes[s]
    write(paste0(1, " ", 1, " ", st, " ", st), fileout, append = TRUE, ...)
    write(paste0(ar[1], " ", ar[2]), fileout, append = TRUE, ...)
    x <- cbind(1, S)
    write(t(x), fileout, ncolumns = ncol(x), append = TRUE, sep = "  \t", ...)
  }
  invisible(NULL)
}

##' Write all data files from a samjr data object
##'
##' Writes the per-(year, age) matrices and the catch / survey
##' observations of a \code{sam_data} list to ICES \code{.dat} files in
##' \code{dir}. Mirrors SAM's \code{write.data.files}; samjr only ships
##' one residual catch fleet, so \code{cn.dat} is always a single file.
##'
##' @param data a \code{sam_data} list (or a fitted \code{sam} object).
##' @param dir output directory.
##' @param writeToOne accepted for API compatibility with SAM; ignored.
##' @param ... arguments forwarded to \code{\link{write}}.
##' @return invisibly \code{NULL}.
##' @export
write.data.files <- function(data, dir = ".", writeToOne = TRUE, ...){
  if(inherits(data, "sam")) data <- data$data
  od <- setwd(dir); on.exit(setwd(od))
  write.ices(data$catchMeanWeight, "cw.dat", writeToOne = writeToOne, ...)
  write.ices(data$disMeanWeight,   "dw.dat", writeToOne = writeToOne, ...)
  write.ices(data$landMeanWeight,  "lw.dat", writeToOne = writeToOne, ...)
  write.ices(data$landFrac,        "lf.dat", writeToOne = writeToOne, ...)
  write.ices(data$propMat,         "mo.dat", writeToOne = writeToOne, ...)
  write.ices(data$stockMeanWeight, "sw.dat", writeToOne = writeToOne, ...)
  write.ices(data$propF,           "pf.dat", writeToOne = writeToOne, ...)
  write.ices(data$propM,           "pm.dat", writeToOne = writeToOne, ...)
  write.ices(data$natMor,          "nm.dat", writeToOne = writeToOne, ...)
  catchFleet <- which(data$fleetTypes == 0)[1]
  write.ices(getFleet(data, catchFleet), "cn.dat", writeToOne = writeToOne, ...)
  write.surveys(data, "survey.dat", ...)
  invisible(NULL)
}
