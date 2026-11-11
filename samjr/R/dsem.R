
#' @title Make path matrices
#'
#' @description Constructs path matrices for dynamic structural equation model (DSEM)
#'        using a vector of parameters and specification of the DSEM
#'
#' @param beta_p vector parameters.
#' @param model matrix or data.frame with the following columns, and one row
#'   per one-headed or two-headed arrow in the dynamic structural model:
#' \describe{
#' \item{direction}{whether a path coefficient is one-headed (1) or two-headed (2)}
#' \item{lag}{whether the lag associated with a given coefficient}
#' \item{start}{starting value, used when \code{parameter=0}}
#' \item{parameter}{The parameter number from \code{beta_p} associated with a given path}
#' \item{first}{The variable at the tail of a given path}
#' \item{second}{The variable at the head of a given path}
#' }
#' @param times integer-vector of times to use when defining matrices
#' @param variables character-vector listing variables
#'
#' @importFrom Matrix solve Diagonal sparseMatrix drop0 kronecker bandSparse
#' @importFrom RTMB ADoverload AD
#'
#' @details
#' When \code{length(times)} is \eqn{T} and \code{length(variables)} is \eqn{J},
#' \code{make_matrices} returns matrices of dimension \eqn{TJ \times TJ} representing
#' paths among \eqn{vec(\mathbf{X})} where matrix \eqn{\mathbf{X}} has dimension
#' \eqn{T \times J} and \eqn{vec} stacks columns into a single long vector
#'
#' @return
#' A named list of matrices including:
#' \describe{
#' \item{P_kk}{The matrix of interactions, i.e., one-headed arrows}
#' \item{G_kk}{The matrix of exogenous covariance, i.e., two-headed arrows}
#' }
#' @export
dsem_make_matrices <-
function( beta_p,
          model,
          times,
          variables ){

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  model <- as.data.frame(model)
  model[,'parameter'] = as.integer(model[,'parameter'])

  # Combine fixed, estimated, and mapped parameters into vector
  beta_i = rep(0, nrow(model))
  off = which(model[,'parameter'] == 0)
  if( length(off) > 0 ){
    beta_i[off] = as.numeric(model[off,'start'])
  }
  not_off = which(model[,'parameter'] > 0)
  if( length(not_off) > 0 ){
    beta_i[not_off] = beta_p[model[not_off,'parameter']]
  }

  # Loop through paths
  P_kk = drop0(sparseMatrix( i=1, j=1, x=0, dims=rep(length(variables)*length(times),2) ))   # Make with a zero
  #P_kk = AD(P_kk)
  G_kk = (P_kk)
  for( i in seq_len(nrow(model)) ){
    lag = as.numeric(model[i,'lag'])
    # Time-lag matrix
    L_tt = bandSparse( n = length(times),
                       m = length(times),
                       k = -1 * lag,  # positive lag means below diagonal
                       diagonals = list(rep(1,length(times))) )

    # Interaction matrix
    P_jj = sparseMatrix( i = match(model[i,'second'],variables),
                         j = match(model[i,'first'],variables),
                         x = 1,
                         dims = rep(length(variables),2) )

    # Assemble
    tmp_kk = kronecker(P_jj, L_tt)
    if(abs(as.numeric(model[i,'direction']))==1){
      P_kk = P_kk + beta_i[i] * tmp_kk # AD(tmp_kk)
    }else{
      G_kk = G_kk + beta_i[i] * tmp_kk # AD(tmp_kk)
    }
  }

  # Diagonal component
  I_kk = Diagonal(nrow(P_kk))

  # Assemble
  IminusP_kk = AD(I_kk - P_kk)
  invV_kk = AD(G_kk)
  invV_kk@x = AD(1 / G_kk@x^2)

  out = list(
    "P_kk" = P_kk,
    "G_kk" = G_kk,
    "invV_kk" = invV_kk,
    "IminusP_kk" = IminusP_kk
  )
  return(out)
}

#' @title Make a RAM (Reticular Action Model)
#'
#' @description \code{read_model} converts SEM arrow notation to \code{model} describing SEM parameters
#'
#' @param sem Specification for time-series structural equation model structure
#'        including lagged or simultaneous effects.  See Details section in
#'        \code{dsem::make_dsem_ram} for more description
#' @param times A character vector listing the set of times in order
#' @param variables A character vector listing the set of variables
#' @param covs A character vector listing variables for which to estimate a standard deviation
#' @param quiet Boolean indicating whether to print messages to terminal
#'
#' @details
#' See \code{\link{make_dsem_ram}} for details
read_model <-
function( sem,
          times,
          variables,
          covs = NULL,
          quiet = FALSE ){

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  ####### Error checks
  if( !is.numeric(times) ) stop("`times` must be numeric in `read_model`")

  ####### Define local functions
  #
  add.variances <- function() {
      variables <- need.variance()
      nvars <- length(variables)
      if (nvars == 0)
          return(model)
      message("NOTE: adding ", nvars, " variances to the model")
      paths <- character(nvars)
      par.names <- character(nvars)
      for (i in 1:nvars) {
          paths[i] <- paste(variables[i], "<->", variables[i])
          par.names[i] <- paste("V[", variables[i], "]", sep = "")
      }
      model.2 <- data.frame(
        'path' = c(model$path, paths),
        'lag' = c(model$lag, rep(0,nvars)),
        'name' = c(model$name, par.names),
        'start' = c(model$start, rep(NA, length(paths))) )
      model.2
  }
  need.variance <- function() {
      all.vars <- classify_variables(model)
      exo.vars <- all.vars$exogenous
      end.vars <- all.vars$endogenous
      variables <- logical(0)
      for (i in seq_len(nrow(model))) {
          paths = model[i,1]
          lag = model[i,2]
          vars <- gsub(pattern=" ", replacement="", x=paths)
          vars <- sub("-*>", "->", sub("<-*", "<-", vars))
          vars <- sub("<->|<-", "->", vars)
          vars <- strsplit(vars, "->")[[1]]
          if ((vars[1] != vars[2]) | (lag != 0)) {
              for (a.variable in vars) {
                if (is.na(variables[a.variable]))
                  variables[a.variable] <- TRUE
              }
          }
          else {
              variables[vars[1]] <- FALSE
          }
      }
      if (!exog.variances && length(exo.vars) > 0)
          variables[exo.vars] <- FALSE
      if (!endog.variances && length(end.vars) > 0)
          variables[end.vars] <- FALSE
      names(variables)[variables]
  }

  ####### Step 2 -- Make RAM
  # convert to data frame
  model = scan( text = sem,
                what = list(path = "", lag = 1, par = "", start = 1, dump = ""),
                sep = ",",
                strip.white = TRUE,
                comment.char = "#",
                fill = TRUE,
                quiet = quiet)
  model$path <- gsub("\\t", " ", model$path)
  model$par[model$par == ""] <- NA
  model <- data.frame( "path" = model$path,
                       "lag" = model$lag,
                       "name" = model$par,
                       "start" = model$start )

  # Adding a SD automatically
  if( !is.null(covs) ){
    for (cov in covs) {
      vars <- strsplit(cov, "[ ,]+")[[1]]
      nvar <- length(vars)
      for (i in 1:nvar) {
      for (j in i:nvar) {
        p1 = paste(vars[i], "<->", vars[j])
        p2 = if (i==j) paste("V[", vars[i], "]", sep = "") else paste("C[",vars[i], ",", vars[j], "]", sep = "")
        p3 = NA
        row <- c(p1, 0, p2, p3)
        if( any((row[1]==model[,1]) & (row[2]==model[,2])) ){
          next
        }else{
          model <- rbind(model, row, deparse.level = 0)
        }
      }}
    }
  }

  exog.variances = endog.variances = TRUE
  model = add.variances()

  # Add parameter column
  par.names = model[, 3]
  pars = na.omit(unique(par.names))
  par.nos = apply(outer(pars, par.names, "=="), 2, which)
  par.nos = unlist(sapply( par.nos, FUN=function(x) ifelse(length(x)==0, 0, x) ))
  if (is.null(par.nos)) par.nos = rep(0, nrow(model))
  model = cbind( model, "parameter"=par.nos )

  # Add incidence to model
  model = cbind( model, first=NA, second=NA, direction=NA )
  for( i in seq_len(nrow(model)) ){
    path = parse_path(model[i,1])
    model[i,c('first','second','direction')] = unlist( path[c('first','second','direction')] )
  }

  return(model)
}


#' @title Parse path
#'
#' @description \code{parse_path} is copied from \code{sem::parse.path}
#'
#' @details
#' Copied from package `sem` under licence GPL (>= 2) with permission from John Fox
#'
#' @return Tagged-list defining variables and direction for a specified path coefficient
#'
#' @param path text to parse
#' @export
parse_path <-
function( path ){
  path.1 <- gsub("-", "", gsub(" ", "", path))
  direction <- if(regexpr("<>", path.1) > 0){
    2
  }else if(regexpr("<", path.1) > 0){
    -1
  }else if(regexpr(">", path.1) > 0){
    1
  }else{
    stop(paste("ill-formed path:", path))
  }
  path.1 <- strsplit(path.1, "[<>]")[[1]]
  out = list(first = path.1[1], second = path.1[length(path.1)], direction = direction)
  return(out)
}


#' @title Classify variables path
#'
#' @description \code{classify_variables} is copied from \code{sem:::classifyVariables}
#'
#' @details
#' Copied from package `sem` under licence GPL (>= 2) with permission from John Fox
#'
#' @param model SEM model
#'
#' @return Tagged-list defining exogenous and endogenous variables
#' @export
classify_variables <-
function( model ){

    variables <- logical(0)
    for (paths in model[, 1]) {
        vars <- gsub(pattern=" ", replacement="", x=paths)
        vars <- sub("-*>", "->", sub("<-*", "<-", vars))
        if (grepl("<->", vars)) {
            vars <- strsplit(vars, "<->")[[1]]
            if (is.na(variables[vars[1]]))
                variables[vars[1]] <- FALSE
            if (is.na(variables[vars[2]]))
                variables[vars[2]] <- FALSE
        }
        else if (grepl("->", vars)) {
            vars <- strsplit(vars, "->")[[1]]
            if (is.na(variables[vars[1]]))
                variables[vars[1]] <- FALSE
            variables[vars[2]] <- TRUE
        }
        else if (grepl("<-", vars)) {
            vars <- strsplit(vars, "<-")[[1]]
            if (is.na(variables[vars[2]]))
                variables[vars[2]] <- FALSE
            variables[vars[1]] <- TRUE
        }
        else stop("incorrectly specified model")
    }
    list(endogenous = names(variables[variables]), exogenous = names(variables[!variables]))
}
