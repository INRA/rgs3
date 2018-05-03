## Copyright 2016-2018 Institut National de la Recherche Agronomique (INRA)
##
## This file is part of rgs3.
##
## rgs3 is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## rgs3 is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public
## License along with rgs3.  If not, see
## <http://www.gnu.org/licenses/>.

##' Check data
##'
##' Return TRUE if the inputs are properly formatted and consistent.
##' @param x data frame
##' @param inds vector
##' @param col.id numeric
##' @param col.traits numeric vector
##' @param binary.traits logical vector
##' @return logical
##' @author Timothee Flutre
##' @export
isValidData <- function(x, inds, col.id, col.traits, binary.traits){
  all(is.data.frame(x),
      ncol(x) >= 2, # at least individual identifiers and trait values
      ifelse(! is.null(inds),
             all(is.vector(inds),
                 is.numeric(inds),
                 ! is.null(names(inds)),
                 all(inds == 1:length(inds))),
             TRUE),
      length(col.id) == 1,
      col.id >= 1,
      col.id <= ncol(x),
      is.factor(x[,col.id]),
      ifelse(! is.null(inds),
             all(levels(x[,col.id]) %in% names(inds)),
             TRUE),
      length(col.traits) >= 1,
      length(col.traits) <= ncol(x) - 1,
      all(col.traits >= 1),
      all(col.traits <= ncol(x)),
      all(! is.factor(x[,col.traits])),
      ! col.id %in% col.traits,
      all(is.logical(binary.traits)),
      length(binary.traits) == length(col.traits),
      ifelse(any(binary.traits),
             all(x[,col.traits[which(binary.traits)]] %in% c(0,1,2,NA)),
             TRUE))
}

##' Write data
##'
##' Write the data file for GS3.
##' @param x data frame containing individual identifiers (names of parameter "inds"), trait values (possibly several traits), covariables and cross-classified factors. Missing data coded as "NA" will be replaced by "0" if it is in a 'trait' column and the trait is binary, and by "-9999" otherwise.
##' @param file path to the text file to which x will be written
##' @param inds named vector with values from 1 to I and names the corresponding individual identifiers; this very same vector should also be used with \code{\link{writeGenosForGs3}}
##' @param col.id column index of the individual identifiers in x
##' @param col.traits column indices of traits in x
##' @param binary.traits logical vector indicating if traits are binary or not (should be of the same length as col.traits)
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{writeGenosForGs3}}, \code{\link{writeConfigForGs3}}
##' @export
writeDataForGs3 <- function(x, file, inds, col.id=1, col.traits=2, binary.traits=FALSE){
  stopifnot(isValidData(x, inds, col.id, col.traits, binary.traits),
            is.character(file))

  ## handle factor levels
  x <- droplevels(x)
  for(c in 1:ncol(x))
    if(c != col.id & ! c %in% col.traits & is.factor(x[,c]))
      levels(x[,c]) <- 1:nlevels(x[,c])

  ## handle genotype identifiers
  x[,col.id] <- as.character(x[,col.id])
  x[,col.id] <- inds[match(x[,col.id], names(inds))]

  ## handle missing values
  for(c in 1:ncol(x)){
    if(c == col.id)
      next
    idx <- which(is.na(x[,c]))
    if(length(idx) > 0){
      if(c %in% col.traits){
        if(binary.traits[which(c == col.traits)]){
          x[idx, c] <- 0
        } else
          x[idx, c] <- -9999
      } else
        x[idx, c] <- -9999
    }
  }

  utils::write.table(x=x,
                     file=file,
                     quote=FALSE,
                     sep=" ",
                     row.names=FALSE,
                     col.names=FALSE)
}

##' Check genotypes
##'
##' Return TRUE if the input corresponds to a properly-formatted genotype matrix.
##' @param x matrix with SNP genotypes
##' @param inds vector
##' @return logical
##' @author Timothee Flutre
##' @export
isValidGenos <- function(x, inds){
  all(is.matrix(x),
      all(x %in% c(0,1,2,NA,5)),
      ! is.null(rownames(x)),
      ifelse(! is.null(inds),
             all(is.vector(inds),
                 is.numeric(inds),
                 ! is.null(names(inds)),
                 all(inds == 1:length(inds)),
                 all(names(inds) %in% rownames(x))),
             TRUE))
}

##' Write genotypes
##'
##' Write the genotype file for GS3.
##' @param x matrix with SNP genotypes encoded as allele dose (i.e. 0/1/2 or NA/5), with individuals in rows and SNPs in columns. Row names are compulsory. Missing data coded as "NA" will be replaced by "5". Missing data already coded as "5" won't be modified.
##' @param file path to the text file to which x will be written
##' @param inds named vector with values from 1 to I and names the corresponding individual identifiers; this very same vector should also be used with \code{\link{writeDataForGs3}}
##' @return nothing
##' @author Timothee Flutre
##' @seealso \code{\link{writeDataForGs3}}, \code{\link{writeConfigForGs3}}
##' @export
writeGenosForGs3 <- function(x, file, inds){
  stopifnot(isValidGenos(x, inds))

  ## handle genotype identifiers
  x <- x[names(inds),]
  rownames(x) <- inds[match(rownames(x), names(inds))]

  ## handle missing values
  x[is.na(x)] <- 5

  ## reformat (a bit slow)
  tmp <- cbind(sprintf("% 47s", rownames(x)),
               t(t(apply(x, 1, paste0, collapse=""))))

  utils::write.table(x=tmp,
                     file=file,
                     quote=FALSE,
                     sep=" ",
                     row.names=FALSE,
                     col.names=FALSE)
}

##' Default configuration
##'
##' Return a list corresponding to a default configuration as used in the vignette.
##' @param data.file path to the text file with the data (if not known yet, use NA instead of NULL)
##' @param genos.file path to the text file with the genotypes (if not known yet, use NA instead of NULL)
##' @param ped.file path to the text file with the pedigree (if not used, use NA or "" instead of NULL)
##' @param nb.snps number of SNPs in the SNP genotype file
##' @param method character with the method to use, among BLUP/MCMCBLUP/VCE/PREDICT
##' @param ptl data frame indicating, for each effect, its position in the data file, its type, and its number of levels
##' @param twc numeric vector containing the column index of the trait of interest in the data file, as well as the column index of its weights (0 if no weight)
##' @param rec.id numeric corresponding to the column index of the genotype identifiers in the data file
##' @param use.mix character indicating by \code{"T"} (or \code{"F"}) if the BayesCPi model should be fitted (or not)
##' @param blasso logical indicating if the Bayesian Lasso model should be fitted (requires \code{method="VCE"} and \code{use.mix="F"})
##' @param task.id character containing the task identifier used as prefix for the output files
##' @return list
##' @author Timothee Flutre
##' @seealso \code{\link{writeConfigForGs3}}, \code{\link{isValidConfig}}
##' @export
getDefaultConfig <- function(data.file=NA, genos.file=NA, ped.file="",
                             nb.snps=NA, method="VCE", ptl=NULL, twc=c(NA,0),
                             rec.id=NA, use.mix="F", blasso=FALSE,
                             task.id="GS3"){
  stopifnot(! is.null(data.file),
            ! is.null(genos.file),
            ! is.null(ped.file))
  if(is.na(ped.file))
    ped.file <- ""
  if(! is.na(nb.snps))
    stopifnot(is.numeric(nb.snps),
              length(nb.snps) == 1,
              nb.snps > 0)
  if(! is.null(ptl)){
    stopifnot(is.data.frame(ptl),
              ncol(ptl) == 3,
              all(colnames(ptl) %in% c("position", "type", "nlevels")),
              all(ptl$type %in% c("cross", "cov", "add_SNP", "dom_SNP",
                                  "add_animal", "perm_diagonal")),
              ifelse("add_animal" %in% ptl$type,
                     ped.file != "", TRUE))
  } else
    ptl <- data.frame(position=NA, type=NA, nlevels=NA)
  if(! is.na(rec.id))
    stopifnot(is.numeric(rec.id),
              length(rec.id) == 1,
              rec.id > 0)
  stopifnot(is.character(task.id),
            length(task.id) == 1)
  stopifnot(is.character(use.mix),
            length(use.mix) == 1,
            use.mix %in% c("T", "F"),
            is.logical(blasso),
            ifelse(blasso,
                   all(method == "VCE", use.mix == FALSE),
                   TRUE))

  config <- list(data.file=data.file,
                 genos.file=genos.file,
                 ped.file=ped.file,
                 num.loci=nb.snps,
                 method=method,
                 simul="F",
                 niter=10000,
                 burnin=2000,
                 thin=10,
                 conv.crit="1d-8",
                 correct=1000,
                 vcs.file=paste0(task.id, "_vcs.txt"),
                 sol.file=paste0(task.id, "_sol.txt"),
                 twc=twc,
                 num.eff=nrow(ptl),
                 ptl=ptl,
                 vc=data.frame(var=c("vara", "vard", "varg", "varp", "vare"),
                               exp=c("2.52d-04","1.75d-06","3.56","2.15","0.19"),
                               df=rep("-2", 5),
                               stringsAsFactors=FALSE),
                 rec.id=rec.id,
                 cont="F",
                 mod=rep("T", nrow(ptl)),
                 ap=c(1,10),
                 dp=c(1,1),
                 use.mix=use.mix,
                 blasso=blasso)

  return(config)
}

##' Check configuration
##'
##'
##' @param config list containing the configuration for GS3
##' @return logical
##' @author Timothee Flutre
##' @seealso \code{\link{getDefaultConfig}}, \code{\link{writeConfigForGs3}}
##' @export
isValidConfig <- function(config){
  if(is.list(config)){
    if(all(c("data.file", "genos.file", "ped.file",
             "num.loci", "method", "simul", "niter", "burnin", "thin",
             "conv.crit", "correct", "vcs.file", "sol.file", "twc", "num.eff",
             "ptl", "vc", "rec.id", "cont", "mod", "ap", "dp", "use.mix",
             "blasso") %in% names(config))){
      if(is.na(config$ped.file))
        config$ped.file <- ""
      all(ifelse(! is.na(config$data.file),
                 file.exists(config$data.file),
                 TRUE),
          ifelse(! is.na(config$genos.file),
                 file.exists(config$genos.file),
                 TRUE),
          ifelse(config$ped.file != "",
                 file.exists(config$ped.file),
                 TRUE),
          ! is.na(config$num.loci),
          config$method %in% c("BLUP", "MCMCBLUP", "VCE", "PREDICT"),
          is.vector(config$twc),
          length(config$twc) == 2,
          all(! is.na(config$twc)),
          config$twc[1] > 0,
          config$twc[2] >= 0,
          config$twc[1] != config$twc[2],
          is.data.frame(config$ptl),
          ncol(config$ptl) == 3,
          all(colnames(config$ptl) %in% c("position", "type", "nlevels")),
          all(config$ptl$type %in% c("cross", "cov", "add_SNP", "dom_SNP",
                                     "add_animal", "perm_diagonal")),
          ifelse("add_animal" %in% config$ptl$type,
                 config$ped.file != "", TRUE),
          nrow(config$ptl) == config$num.eff,
          is.vector(config$mod),
          length(config$mod) == config$num.eff,
          all(config$mod %in% c("T", "F")),
          is.numeric(config$burnin),
          length(config$burnin) == 1,
          config$burnin >= 0,
          is.numeric(config$niter),
          length(config$niter) == 1,
          config$niter >= 1,
          is.numeric(config$thin),
          length(config$thin) == 1,
          config$thin >= 0,
          config$thin < config$niter,
          is.data.frame(config$vc),
          ncol(config$vc) == 3,
          all(colnames(config$vc) %in% c("var", "exp", "df")),
          ! is.na(config$rec.id),
          is.vector(config$ap),
          length(config$ap) == 2,
          is.vector(config$dp),
          length(config$dp) == 2,
          is.character(config$use.mix),
          length(config$use.mix) == 1,
          config$use.mix %in% c("T", "F"),
          is.logical(config$blasso),
          length(config$blasso) == 1,
          ifelse(config$blasso,
                 all(config$method == "VCE", config$use.mix == "F"),
                 TRUE)
          )
    }
  } else
    FALSE
}

##' Write configuration
##'
##' Write the configuration file for GS3.
##' @param config list containing the configuration for GS3 via the following components:
##' \describe{
##'   \item{num.loci}{number of loci}
##'   \item{method}{BLUP/MCMCBLUP/VCE/PREDICT}
##'   \item{simul}{one-letter character indicating if simulations should be performed (T) or not (F)}
##'   \item{niter}{number of iterations for the Gibbs sampler}
##'   \item{burnin}{burn-in for the Gibbs sampler}
##'   \item{thin}{thinning for the Gibbs sampler}
##'   \item{conv.crit}{convergence criterion (meaningful if BLUP)}
##'   \item{correct}{correction (to avoid numerical problems)}
##'   \item{vcs.file}{path to the text file to which the variance component samples will be written}
##'   \item{sol.file}{path to the text file to which the solutions will be written}
##'   \item{twc}{2-element vector which first element corresponds to the column index of the trait values in the data file, and the second to the column index of the weights in the data file (use 0 if no weight)}
##'   \item{num.eff}{number of effects}
##'   \item{ptl}{3-column data frame indicating, for each covariable/factor, its position (column) in the data file, type of effect and number of levels}
##'   \item{vc}{3-column data frame indicating, for each variance component, its expected value and degrees of freedom}
##'   \item{rec.id}{1-element vector indicating the column in the data file corresponding to the genotype identifiers}
##'   \item{cont}{vector with T or F indicating if the MCMC run is a continuation of a previous, interrupted one}
##'   \item{mod}{vector with T or F for each covariable/factor indicating if it has to be included or not in the model}
##'   \item{ap}{prior proportions of the BayesCPi mixture}
##'   \item{dp}{prior proportions of the BayesCPi mixture}
##'   \item{use.mix}{one-letter character indicating if the Bayes C pi prior should be used (T) or not (F)}
##'   \item{blasso}{boolean indicating if the Bayesian lasso prior should be used or not}
##' }
##' @param task.id character containing the task identifier used as prefix for the configuration file; will be followed by \code{"_config.txt"}
##' @return path to the text file to which the configuration for GS3 will be written
##' @author Timothee Flutre
##' @seealso \code{\link{getDefaultConfig}}, \code{\link{execGs3}}
##' @export
writeConfigForGs3 <- function(config,
                              task.id="GS3"){
  stopifnot(isValidConfig(config),
            file.exists(config$data.file))
  tmp <- utils::read.table(file=config$data.file, sep=" ", nrows=2)
  stopifnot(config$twc[1] <= ncol(tmp),
            ifelse(config$twc[2] != 0,
                   config$twc[2] <= ncol(tmp),
                   TRUE))
  if(config$ped.file != "")
    stopifnot(file.exists(config$ped.file))
  if(! is.na(config$genos.file)){
    stopifnot(file.exists(config$genos.file))
    tmp <- readLines(con=config$genos.file, n=1)
    stopifnot(config$num.loci == nchar(tmp) - 48) # see writeGenosForGs3()
  }

  config.file <- paste0(task.id, "_config.txt")

  txt <- paste0("DATAFILE",
                "\n", config$data.file)
  txt <- paste0(txt, "\nPEDIGREE FILE",
                "\n", config$ped.file)
  txt <- paste0(txt, "\nGENOTYPE FILE",
                "\n", config$genos.file)
  txt <- paste0(txt, "\nNUMBER OF LOCI (might be 0)",
                "\n", sprintf("%i", config$num.loci))
  txt <- paste0(txt, "\nMETHOD (BLUP/MCMCBLUP/VCE/PREDICT)",
                "\n", config$method)
  txt <- paste0(txt, "\nSIMULATION",
                "\n", config$simul)
  txt <- paste0(txt, "\nGIBBS SAMPLING PARAMETERS")
  txt <- paste0(txt, "\nNITER",
                "\n", sprintf("%i", config$niter))
  txt <- paste0(txt, "\nBURNIN",
                "\n", sprintf("%i", config$burnin))
  txt <- paste0(txt, "\nTHIN",
                "\n", sprintf("%i", config$thin))
  txt <- paste0(txt, "\nCONV_CRIT (MEANINGFUL IF BLUP)",
                "\n", config$conv.crit)
  txt <- paste0(txt, "\nCORRECTION (to avoid numerical problems)",
                "\n", config$correct)
  txt <- paste0(txt, "\nVARIANCE COMPONENTS SAMPLES",
                "\n", config$vcs.file)
  txt <- paste0(txt, "\nSOLUTION FILE",
                "\n", config$sol.file)
  txt <- paste0(txt, "\nTRAIT AND WEIGHT COLUMNS",
                "\n", paste0(config$twc, collapse=" "))
  txt <- paste0(txt, "\nNUMBER OF EFFECTS",
                "\n", config$num.eff)
  txt <- paste0(txt, "\nPOSITION IN DATA FILE TYPE OF EFFECT  NUMBER OF LEVELS")
  for(i in 1:nrow(config$ptl))
    txt <- paste0(txt, "\n", config$ptl$position[i], " ", config$ptl$type[i],
                  " ", config$ptl$nlevels[i])
  txt <- paste0(txt, "\nVARIANCE COMPONENTS (fixed for any BLUP, starting values for VCE)")
  for(i in 1:nrow(config$vc))
    txt <- paste0(txt, "\n", config$vc$var[i],
                  "\n", config$vc$exp[i], " ", config$vc$df[i])
  txt <- paste0(txt, "\nRECORD ID",
                "\n", config$rec.id)
  txt <- paste0(txt, "\nCONTINUATION (T/F)",
                "\n", config$cont)
  txt <- paste0(txt, "\nMODEL (T/F for each effect)",
                "\n", paste0(config$mod, collapse=" "))
  txt <- paste0(txt, "\nA PRIORI a",
                "\n", paste0(config$ap, collapse=" "))
  txt <- paste0(txt, "\na PRIORI D",
                "\n", paste0(config$dp, collapse=" "))
  txt <- paste0(txt, "\nUSE MIXTURE",
                "\n", config$use.mix)
  if(config$blasso)
    txt <- paste0(txt, "\n#OPTION BayesianLasso Tibshirani")
  txt <- paste0(txt, "\n")

  cat(txt, file=config.file)

  return(config.file)
}

##' Execute GS3
##'
##' Execute GS3 via a system call.
##' The expected executable name in the PATH is \code{gs3.exe} on Windows, \code{gs3} otherwise.
##' @param config.file path to the text file containing the configuration for GS3
##' @param task.id character containing the task identifier used as prefix for the "stdout/stderr" and "freq" output files
##' @return path to the text file to which the stdout and stderr will be written
##' @author Timothee Flutre
##' @seealso \code{\link{writeConfigForGs3}}, \code{\link{vcs2mcmc}}, \code{\link{cleanGs3}}
##' @export
execGs3 <- function(config.file, task.id="GS3"){
  stopifnot(file.exists(config.file),
            is.character(task.id))

  stdouterr.file <- paste0(task.id, "_stdouterr.txt")

  if(Sys.info()["sysname"] == "Windows"){
    cmd <- "gs3.exe"
    ret <- system2(command=cmd, args=c(config.file),
                   stdout=stdouterr.file, stderr=stdouterr.file, wait=TRUE,
                   invisible=FALSE)
  } else{
    cmd <- "gs3"
    ret <- system2(command=cmd, args=c(config.file),
                   stdout=stdouterr.file, stderr=stdouterr.file, wait=TRUE)
  }

  if(ret != 0){
    msg <- paste0("'", cmd, "' returned '", ret, "'")
    warning(msg)
  }

  return(stdouterr.file)
}

identifyOptModComps <- function(config){
  stopifnot(is.list(config),
            all(c("ptl","mod","ped.file","use.mix","blasso") %in%
                names(config)))

  out <- rep(FALSE, 5)
  names(out) <- c("has.d", "has.g", "has.p", "is.bayesCpi", "is.blasso")

  idx <- which(config$ptl$type == "dom_SNP")
  if(length(idx) == 1)
    out["has.d"] <- config$mod[idx] == "T"

  if(config$ped.file != ""){
    idx <- which(config$ptl$type == "add_animal")
    if(length(idx) == 1)
      out["has.g"] <- config$mod[idx] == "T"
  }

  idx <- which(config$ptl$type == "perm_diagonal")
  if(length(idx) == 1)
    out["has.p"] <- config$mod[idx] == "T"

  out["is.bayesCpi"] <- config$use.mix == "T"

  out["is.blasso"] <- config$blasso

  return(out)
}

addStatGenoVarComp <- function(dat, afs, has.d=FALSE, is.bayesCpi=FALSE){
  stopifnot(is.data.frame(dat),
            "vara" %in% colnames(dat),
            is.numeric(afs),
            all(! is.na(afs)),
            all(afs >= 0),
            all(afs <= 1),
            is.logical(has.d),
            is.logical(is.bayesCpi))
  if(has.d)
    stopifnot("vard" %in% colnames(dat))
  if(is.bayesCpi){
    stopifnot("pa_1" %in% colnames(dat))
    if(has.d)
      stopifnot("pd_1" %in% colnames(dat))
  }

  out <- dat

  ## see Vitezica et al (2013): page 1226, bottom of the right column
  ## as well as the GS3 manual (section "Use", subsection "output")
  out$varA <- out$vara * 2 * sum(afs * (1 - afs))
  if(is.bayesCpi)
    out$varA <- out$varA * out$pa_1

  if(has.d){
    tmp <- out$vard * 2 * sum(afs * (1 - afs) * (1 - 2 * afs)^2)
    if(is.bayesCpi)
      tmp <- tmp * out$pd_1
    out$varA <- out$varA + tmp

    out$varD <- out$vard * 2^2 * sum(afs^2 * (1 - afs)^2)
    if(is.bayesCpi)
      out$varD <- out$varD * out$pd_1
  }

  return(out)
}

##' Load GS3 results
##'
##' Read the file containing the variance components' samples into a \code{\link[coda]{mcmc.list}} object.
##' @param config list containing the configuration for GS3
##' @param afs vector of allele frequencies (of all the SNPs which genotypes were provided, and only them); used to compute the variance of additive (and dominance) genotypic values from the variance of additive (and dominance) SNP effects (see Vitezica et al, 2013)
##' @return \code{\link[coda]{mcmc.list}}
##' @author Timothee Flutre
##' @seealso \code{\link{execGs3}}
##' @examples
##' \dontrun{vcs <- vcs2mcmc(config$vcs.file)
##' summary(vcs)
##' coda::effectiveSize(vcs)
##' genos <- as.matrix(read.table(genos.file))
##' afs <- colMeans(genos) / 2
##' vcs <- vcs2mcmc(config$vcs.file, afs)
##' summary(vcs)}
##' @export
vcs2mcmc <- function(config, afs=NULL){
  if(! requireNamespace("coda", quietly=TRUE))
    stop("Pkg coda needed for this function to work. Please install it.",
         call.=FALSE)
  stopifnot(isValidConfig(config))
  if(! is.null(afs))
    stopifnot(is.numeric(afs),
              all(! is.na(afs)),
              all(afs >= 0),
              all(afs <= 1))

  ## load the samples
  d <- utils::read.table(config$vcs.file, header=TRUE, check.names=FALSE)

  ## identify the optional model components
  opt.mod.comps <- identifyOptModComps(config)

  ## discard useless columns
  if(! opt.mod.comps["has.d"] & "vard" %in% colnames(d))
    d$vard <- NULL
  if(! opt.mod.comps["has.g"] & "varg" %in% colnames(d))
    d$varg <- NULL
  if(! opt.mod.comps["has.p"] & "varp" %in% colnames(d))
    d$varp <- NULL
  if(! opt.mod.comps["is.bayesCpi"]){
    if("pa_1" %in% colnames(d))
      d$pa_1 <- NULL
    if("pd_1" %in% colnames(d))
      d$pd_1 <- NULL
  } else{
    if(! opt.mod.comps["has.d"] & "pd_1" %in% colnames(d))
      d$pd_1 <- NULL
  }
  if(! opt.mod.comps["is.blasso"] & "lambda2" %in% colnames(d))
    d$lambda2 <- NULL
  for(j in seq_along(d))
    if(! is.numeric(d[[j]]))
      d[[j]] <- NULL

  ## compute the variances of add and dom genotypic values
  if(! is.null(afs))
    d <- addStatGenoVarComp(dat=d, afs=afs,
                            has.d=opt.mod.comps["has.d"],
                            is.bayesCpi=opt.mod.comps["is.bayesCpi"])

  ## compute narrow-sense heritability (h^2)
  d$h2 <- d$varA / (d$varA +
                    ifelse(! is.null(d$varD), d$varD, 0) +
                    ifelse(! is.null(d$varg), d$varg, 0) +
                    ifelse(! is.null(d$varp), d$varp, 0) +
                    d$vare)

  ## convert to 'coda' format
  vcs <- coda::mcmc.list(coda::mcmc(d,
                                    start=config$burnin + 1,
                                    end=config$burnin + config$niter,
                                    thin=config$thin))

  return(vcs)
}

##' Clean output files
##'
##' Remove the various output files generated by GS3.
##' @param config list containing the configuration of the GS3 run to be cleaned
##' @param config.file path to the text file to which the configuration for GS3 (corresponding to \code{config}) is written
##' @param task.id character containing the task identifier used as prefix for the output files
##' @author Timothee Flutre
##' @export
cleanGs3 <- function(config, config.file, task.id="GS3"){
  stopifnot(is.list(config),
            file.exists(config.file),
            is.character(task.id))

  ## build the list of files to remove
  files.tormv <- c(config$vcs.file,
                   config$sol.file,
                   config.file,
                   paste0(config.file, "_freq"),
                   paste0(config.file, "_EBVs"),
                   paste0(config.file, "_cont"),
                   paste0(config.file, "_finalEstimates"))
  files.tormv <- c(files.tormv, paste0(task.id, "_stdouterr.txt"))

  ## remove the files
  for(f in files.tormv)
    if(file.exists(f))
      file.remove(f)
}

##' Partition genotypes for cross-validation
##'
##' Return a vector which content corresponds to fold indices and names to genotypes.
##' The fold index of a given genotype indicates that, for this fold, the genotype won't be used for training but for validation.
##' @param geno.names vector of genotype names
##' @param nb.folds number of folds
##' @param seed if not NULL, this seed for the pseudo-random number generator will be used to shuffle genotypes before partitioning per fold
##' @return vector
##' @author Timothee Flutre
##' @seealso \code{\link{crossValWithGs3}}
##' @export
getPartitionGenos <- function(geno.names, nb.folds=10, seed=NULL){
  stopifnot(is.vector(geno.names),
            nb.folds <= length(geno.names),
            ifelse(! is.null(seed), is.numeric(seed), TRUE))

  if(! is.null(seed)){
    set.seed(seed)
    geno.names <- sample(x=geno.names)
  }

  nb.genos <- length(geno.names)
  nb.valid.genos.fold1 <- ceiling(nb.genos / nb.folds)
  fold.indices <- seq(from=1, to=nb.genos, by=nb.valid.genos.fold1)
  valid.geno.names.per.fold <- lapply(fold.indices, function(i){
    geno.names[i:min(i + nb.valid.genos.fold1 - 1, nb.genos)]
  })

  nb.valid.genos.per.fold <- sapply(valid.geno.names.per.fold, length)
  valid.geno.idx.per.fold <- stats::setNames(
      object=rep.int(1:nb.folds, times=nb.valid.genos.per.fold),
      nm=unlist(valid.geno.names.per.fold))

  return(valid.geno.idx.per.fold)
}

##' Single fold of a cross-validation
##'
##' Perform a single fold of a K-fold cross-validation with GS3.
##' @param task.id character containing the task identifier used as prefix for the output files (the fold index will be added)
##' @param rep.id identifier of the current replicate
##' @param fold.ids vector of fold identifiers
##' @param fold.id identifier of the current fold
##' @param valid.geno.idx.per.fold vector indicating to which fold a given genotypes belongs (see \code{\link{getPartitionGenos}})
##' @param dat data frame with phenotypes
##' @param col.id column identifier of the records in \code{dat}
##' @param col.trait column identifier of the trait in \code{dat}
##' @param binary.trait logical indicating if the trait is binary or not
##' @param genos matrix of SNP genotypes
##' @param config list containing the generic configuration for GS3
##' @param ped.file path to the file containing the pedigree
##' @param afs vector of allele frequencies which names are SNP identifiers (column names of \code{genos}); used to compute the variance of additive genotypic values from the variance of additive SNP effects (see Vitezica et al, 2013), then used to compute narrow-sense heritability
##' @param remove.files remove files per fold (none/some/all); use \code{"some"} in real-life applications in order to keep estimates of SNP effects per fold, thereby allowing to perform genomic prediction afterwards by averaging them
##' @param verbose verbosity level (0/1/2)
##' @return named vector with metrics
##' @author Timothee Flutre
##' @seealso \code{\link{crossValWithGs3}}
##' @export
crossValFold <- function(task.id, rep.id, fold.ids, fold.id,
                         valid.geno.idx.per.fold,
                         dat, col.id, col.trait, binary.trait,
                         genos, config, ped.file, afs,
                         remove.files, verbose=1){
  ## identify the optional random variables
  opt.mod.comps <- identifyOptModComps(config)

  ## prepare the output
  out <- stats::setNames(c(fold.id, rep(NA, 4)),
                         c("fold", "genos.train", "obs.train",
                           "var.a.mean", "var.a.sd"))
  if(opt.mod.comps["has.d"]){
    out["var.d.mean"] <- NA
    out["var.d.sd"] <- NA
  }
  out["var.A.mean"] <- NA
  out["var.A.sd"] <- NA
  if(opt.mod.comps["has.d"]){
    out["var.D.mean"] <- NA
    out["var.D.sd"] <- NA
  }
  if(opt.mod.comps["has.g"]){
    out["var.g.mean"] <- NA
    out["var.g.sd"] <- NA
  }
  if(opt.mod.comps["has.p"]){
    out["var.p.mean"] <- NA
    out["var.p.sd"] <- NA
  }
  out["var.e.mean"] <- NA
  out["var.e.sd"] <- NA
  out["h2.mean"] <- NA
  out["h2.sd"] <- NA
  out["genos.valid"] <- NA
  out["obs.valid"] <- NA
  out["rmspe"] <- NA
  out["cor.p"] <- NA
  out["cor.p.div.h"] <- NA
  out["cor.s"] <- NA
  out["reg.intercept"] <- NA
  out["reg.slope"] <- NA

  nb.folds <- length(fold.ids)

  ## prepare temporary file names
  tif <- paste0(task.id, "_r", rep.id, "_f", fold.ids[fold.id])
  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ", fold ", fold.id, "/", nb.folds,
                  ": set up file names...\n", tif)
    write(msg, stdout())
  }
  dat.train.file <- paste0(tif, "_dat-train.txt")
  geno.train.file <- paste0(tif, "_geno-train.txt")
  dat.valid.file <- paste0(tif, "_dat-valid.txt")
  geno.valid.file <- paste0(tif, "_geno-valid.txt")

  ## -------------------------------------------------------------------------
  ## training
  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ", fold ", fold.id, "/", nb.folds,
                  ": start training")
    write(msg, stdout())
  }
  train.genos <- names(which(valid.geno.idx.per.fold != fold.id))
  out["genos.train"] <- length(train.genos)
  dat.train <- droplevels(dat[dat[,col.id] %in% train.genos,])
  out["obs.train"] <- sum(! is.na(dat.train[, col.trait]))
  genos.train <- genos[train.genos,]
  inds.train <- stats::setNames(object=1:nlevels(dat.train[, col.id]),
                                nm=levels(dat.train[, col.id]))
  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ", fold ", fold.id, "/", nb.folds, ":",
                  " nb.genos=", nrow(genos.train),
                  " nb.snps=", ncol(genos.train),
                  " nb.obs=", sum(! is.na(dat.train[, col.trait])))
    write(msg, stdout())
  }
  writeDataForGs3(x = dat.train,
                  file = dat.train.file,
                  inds = inds.train,
                  col.id = col.id,
                  col.traits = col.trait,
                  binary.traits = binary.trait)
  writeGenosForGs3(x = genos.train,
                   file = geno.train.file,
                   inds = inds.train)
  config.train <- config
  config.train$data.file <- dat.train.file
  config.train$genos.file <- geno.train.file
  config.train$ped.file <- ped.file
  config.train$vcs.file <- sub(task.id, tif, config$vcs.file)
  config.train$sol.file <- sub(task.id, tif, config$sol.file)
  config.train.file <- writeConfigForGs3(config = config.train,
                                         task.id = tif)
  stdouterr.train.file <- execGs3(config.file = config.train.file,
                                  task.id = tif)
  if(FALSE){ # for debugging purposes
    stdouterr.train <- readLines(stdouterr.train.file)
    sols <- utils::read.table(file=config.train$sol.file, header=TRUE)
  }
  vcs <- utils::read.table(file=config.train$vcs.file, header=TRUE,
                           check.names=FALSE)
  ## see Vitezica et al (2013)
  out["var.a.mean"] <- mean(vcs$vara)
  out["var.a.sd"] <- stats::sd(vcs$vara)
  if(opt.mod.comps["has.d"]){
    out["var.d.mean"] <- mean(vcs$vard)
    out["var.d.sd"] <- stats::sd(vcs$vard)
  }
  vcs <- addStatGenoVarComp(dat=vcs, afs=afs,
                            has.d=opt.mod.comps["has.d"],
                            is.bayesCpi=opt.mod.comps["is.bayesCpi"])
  out["var.A.mean"] <- mean(vcs$varA)
  out["var.A.sd"] <- stats::sd(vcs$varA)
  if(opt.mod.comps["has.d"]){
    out["var.D.mean"] <- mean(vcs$varD)
    out["var.D.sd"] <- stats::sd(vcs$varD)
  }
  if(opt.mod.comps["has.g"]){
    out["var.g.mean"] <- mean(vcs$varg)
    out["var.g.sd"] <- stats::sd(vcs$varg)
  }
  if(opt.mod.comps["has.p"]){
    out["var.p.mean"] <- mean(vcs$varp)
    out["var.p.sd"] <- stats::sd(vcs$varp)
  }
  out["var.e.mean"] <- mean(vcs$vare)
  out["var.e.sd"] <- stats::sd(vcs$vare)
  vcs$h2 <- vcs$varA / (vcs$varA +
                        ifelse(opt.mod.comps["has.d"], vcs$varD, 0) +
                        ifelse(opt.mod.comps["has.g"], vcs$varg, 0) +
                        ifelse(opt.mod.comps["has.p"], vcs$varp, 0) +
                        vcs$vare)
  out["h2.mean"] <- mean(vcs[, "h2"])
  out["h2.sd"] <- stats::sd(vcs[, "h2"])
  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ", fold ", fold.id, "/", nb.folds,
                  ": end training")
    write(msg, stdout())
  }

  ## -------------------------------------------------------------------------
  ## validation
  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ", fold ", fold.id, "/", nb.folds,
                  ": start validation")
    write(msg, stdout())
  }
  valid.genos <- names(which(valid.geno.idx.per.fold == fold.id))
  stopifnot(all(! valid.genos %in% train.genos))
  out["genos.valid"] <- length(valid.genos)
  dat.valid <- droplevels(dat[dat[,col.id] %in% valid.genos,])
  out["obs.valid"] <- nrow(dat.valid)
  genos.valid <- genos[valid.genos,]
  inds.valid <- stats::setNames(object=1:nlevels(dat.valid[, col.id]),
                                nm=levels(dat.valid[, col.id]))
  dat.valid.to.predict <- dat.valid
  dat.valid.to.predict[, col.trait] <- NA
  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ", fold ", fold.id, "/", nb.folds, ":",
                  " nb.genos=", nrow(genos.valid),
                  " nb.snps=", ncol(genos.valid),
                  " nb.pred=", nrow(dat.valid))
    write(msg, stdout())
  }
  writeDataForGs3(x = dat.valid.to.predict,
                  file = dat.valid.file,
                  inds = inds.valid,
                  col.id = col.id,
                  col.traits = col.trait,
                  binary.traits = binary.trait)
  writeGenosForGs3(x = genos.valid,
                   file = geno.valid.file,
                   inds = inds.valid)
  config.valid <- config.train
  config.valid$data.file <- dat.valid.file
  config.valid$genos.file <- geno.valid.file
  config.valid$ped.file <- ped.file
  config.valid$method <- "PREDICT"
  config.valid.file <- writeConfigForGs3(config = config.valid,
                                         task.id = tif)
  stdouterr.valid.file <- execGs3(config.file = config.valid.file,
                                  task.id = tif)
  if(FALSE) # for debugging purposes
    stdouterr.valid <- readLines(stdouterr.valid.file)
  pred.file <- paste0(config.valid.file, "_predictions")
  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ", fold ", fold.id, "/", nb.folds,
                  ": end validation")
    write(msg, stdout())
  }

  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ", fold ", fold.id, "/", nb.folds,
                  ": assess accuracy...")
    write(msg, stdout())
  }
  if(FALSE){ # for debugging purposes
    ebvs <- utils::read.table(file=paste0(config.valid.file, "_EBVs"),
                              header=TRUE)
  }
  pred <- utils::read.table(file=pred.file, header=TRUE)
  stopifnot(nrow(pred) == nrow(dat.valid))
  stopifnot(colnames(pred) == c("id", "true", "prediction"))
  pred$geno <- names(inds.valid)[match(pred$id, inds.valid)]
  pred <- pred[order(pred$geno),]
  dat.valid <- dat.valid[order(dat.valid[,col.id]),]
  stopifnot(all(pred$geno == as.character(dat.valid[,col.id])))
  ## note that some phenotypes (dat.valid[, col.trait]) can be missing
  out["rmspe"] <- sqrt(mean((dat.valid[, col.trait] -
                             pred$prediction)^2, na.rm=TRUE))
  out["cor.p"] <- stats::cor(x=dat.valid[, col.trait],
                             y=pred$prediction,
                             use="complete.obs",
                             method="pearson")
  out["cor.p.div.h"] <- stats::cor(x=dat.valid[, col.trait],
                                   y=pred$prediction,
                                   use="complete.obs",
                                   method="pearson") /
    sqrt(out["h2.mean"])
  out["cor.s"] <- stats::cor(x=dat.valid[, col.trait],
                             y=pred$prediction,
                             use="complete.obs",
                             method="spearman")
  fit <- stats::lm(dat.valid[, col.trait] ~ pred$prediction)
  out["reg.intercept"] <- stats::coefficients(fit)[1]
  out["reg.slope"] <- stats::coefficients(fit)[2]
  if(FALSE){ # for debugging purposes
    plot(x=pred$prediction, y=dat.valid[, col.trait],
         main=paste0("fold ", fold.id))
    abline(a=0, b=1, lty=2)
    abline(v=mean(pred$prediction), lty=2)
    abline(h=mean(dat.valid[, col.trait]), lty=2)
    abline(fit, col="red")
  }

  ## remove temporary files
  if(remove.files != "none"){
    if(remove.files %in% c("some", "all")){
      files.tormv <- c(dat.train.file,
                       geno.train.file,
                       config.train.file,
                       stdouterr.train.file,
                       config.train$vcs.file,
                       paste0(config.train.file, "_freq"),
                       paste0(config.train.file, "_EBVs"),
                       paste0(config.train.file, "_cont"),
                       dat.valid.file,
                       geno.valid.file,
                       config.valid.file,
                       paste0(config.valid.file, "_freq"),
                       paste0(config.valid.file, "_EBVs"),
                       paste0(config.valid.file, "_cont"),
                       paste0(config.valid.file, "_predictions"),
                       stdouterr.valid.file)
    }
    if(remove.files == "all"){
      files.tormv <- c(files.tormv,
                       paste0(config.train.file, "_finalEstimates"),
                       config.train$sol.file)
    }
    for(f in files.tormv)
      if(file.exists(f))
        file.remove(f)
  }

  return(out)
}

##' K-fold cross-validation
##'
##' Perform K-fold cross-validation with GS3, and report metrics as advised in \href{http://dx.doi.org/10.1534/genetics.112.147983}{Daetwyler et al. (2013)}.
##' Files are saved in the current directory.
##' @param genos matrix of SNP genotypes
##' @param dat data frame with phenotypes
##' @param config list containing the configuration for GS3
##' @param task.id character containing the task identifier used as prefix for the output files (for each fold, its index will be added)
##' @param binary.trait logical
##' @param ped.file path to the file containing the pedigree (if not used, use NA or "" instead of NULL)
##' @param afs vector of allele frequencies which names are SNP identifiers (column names of \code{genos}); if NULL, will be estimated from \code{genos}; used to compute the variance of additive genotypic values from the variance of additive SNP effects (see Vitezica et al, 2013), then used to compute narrow-sense heritability
##' @param rep.id identifier of the current replicate
##' @param nb.folds number of folds
##' @param seed if not NULL, this seed for the pseudo-random number generator will be used to shuffle genotypes before partitioning per fold
##' @param remove.files remove files per fold (none/some/all); use \code{"some"} in real-life applications in order to keep estimates of SNP effects per fold, thereby allowing to perform genomic prediction afterwards by averaging them
##' @param nb.cores number of cores to launch folds in parallel (via \code{\link[parallel]{mclapply}}, on Unix-like computers, or \code{\link[parallel]{parLapply}} on Windows); you can also use \code{\link[parallel]{detectCores}}
##' @param cl object returned by \code{\link[parallel]{makeCluster}}, necessary only if \code{nb.cores > 1} and the computer runs Windows; if NULL in such cases, will be created silently
##' @param verbose verbosity level (0/1/2); there will be a progress bar only for \code{verbose=1}
##' @return data frame
##' @author Helene Muranty, Timothee Flutre
##' @seealso \code{\link{crossValRepWithGs3}}, \code{\link{getDefaultConfig}}, \code{\link{writeConfigForGs3}}, \code{\link{execGs3}}, \code{\link{getPartitionGenos}}
##' @export
crossValWithGs3 <- function(genos,
                            dat,
                            config,
                            task.id="GS3",
                            binary.trait=FALSE,
                            ped.file="",
                            afs=NULL,
                            rep.id=1,
                            nb.folds=10,
                            seed=NULL,
                            remove.files="some",
                            nb.cores=1,
                            cl=NULL,
                            verbose=1){
  requireNamespace("parallel")
  stopifnot(is.list(config),
            all(c("rec.id", "twc", "num.loci") %in%names(config)))
  col.id <- config$rec.id
  col.trait <- config$twc[1]
  stopifnot(! is.null(ped.file))
  if(is.na(ped.file))
    ped.file <- ""
  stopifnot(isValidGenos(genos, NULL),
            config$num.loci == ncol(genos),
            isValidData(dat, NULL, col.id, col.trait, binary.trait),
            all(rownames(genos) %in% levels(dat[, col.id])),
            all(levels(dat[, col.id]) %in% rownames(genos)),
            nb.folds <= nrow(genos),
            remove.files %in% c("none", "some", "all"),
            is.character(task.id),
            length(task.id) == 1,
            is.numeric(nb.cores),
            nb.cores > 0)

  wasClCreated <- FALSE
  if(all(nb.cores > 1, Sys.info()["sysname"] == "Windows", is.null(cl))){
    wasClCreated <- TRUE
    cl <- parallel::makeCluster(nb.cores)
  }

  if(is.null(afs)){
    if(verbose > 0){
      msg <- paste0("rep ", rep.id, ": estimate allele frequencies...")
      write(msg, stdout())
    }
    afs <- colMeans(genos, na.rm=TRUE) / 2
    if(verbose > 0)
      print(summary(afs))
  } else{
    stopifnot(is.numeric(afs),
              ! is.null(names(afs)),
              all(colnames(genos) %in% names(afs)))
    afs <- afs[colnames(genos)]
    stopifnot(all(! is.na(afs)),
              all(afs >= 0),
              all(afs <= 1))
  }

  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ": prepare partitions...")
    write(msg, stdout())
  }
  valid.geno.idx.per.fold <- getPartitionGenos(geno.names=rownames(genos),
                                               nb.folds=nb.folds,
                                               seed=seed)
  stopifnot(! any(duplicated(names(valid.geno.idx.per.fold))))

  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ": perform cross-validation...")
    write(msg, stdout())
  }
  fold.ids <- sprintf(fmt=paste0("%0", floor(log10(nb.folds))+1, "i"),
                      1:nb.folds)
  if(nb.cores == 1){
    out <- lapply(seq(nb.folds), function(fold.id){
      crossValFold(task.id, rep.id, fold.ids, fold.id,
                   valid.geno.idx.per.fold,
                   dat, col.id, col.trait, binary.trait,
                   genos, config, ped.file, afs,
                   remove.files, verbose - 1)
    })
  } else{ # nb.cores > 1
    if(Sys.info()["sysname"] == "Windows"){
      out <- parallel::parLapply(cl, seq(nb.folds), function(fold.id){
        crossValFold(task.id, rep.id, fold.ids, fold.id,
                     valid.geno.idx.per.fold,
                     dat, col.id, col.trait, binary.trait,
                     genos, config, ped.file, afs,
                     remove.files, verbose - 1)
      })
    } else{
      out <- parallel::mclapply(seq(nb.folds), function(fold.id){
        tryCatch({
          crossValFold(task.id, rep.id, fold.ids, fold.id,
                       valid.geno.idx.per.fold,
                       dat, col.id, col.trait, binary.trait,
                       genos, config, ped.file, afs,
                       remove.files, verbose - 1)
        },
        error=function(e) {
          print(e)
          stop(e)
        })
      }, mc.cores=nb.cores)
    } # end of "non-Windows" case
  } # end of "nb.cores > 1" case

  if(verbose > 0){
    msg <- paste0("rep ", rep.id, ": prepare output...")
    write(msg, stdout())
  }
  out <- do.call(rbind, out)
  if(! is.null(seed))
    attr(out, "seed") <- seed

  if(wasClCreated)
    parallel::stopCluster(cl)

  return(out)
}

##' Replicated K-fold cross-validation
##'
##' Perform replicated K-fold cross-validation with GS3, i.e. call \code{\link{crossValWithGs3}} several times, with different seeds.
##' @param genos matrix of SNP genotypes
##' @param dat data frame with phenotypes
##' @param config list containing the configuration for GS3
##' @param task.id character containing the task identifier used as prefix for the output files (for each fold, its index will be added)
##' @param binary.trait logical
##' @param ped.file path to the file containing the pedigree (if not used, use NA or "" instead of NULL)
##' @param afs vector of allele frequencies which names are SNP identifiers (column names of \code{genos}); if NULL, will be estimated from \code{genos}; used to compute the variance of additive genotypic values from the variance of additive SNP effects (see Vitezica et al, 2013), then used to compute narrow-sense heritability
##' @param nb.reps number of replicates (a set of folds will be sampled for each replicate)
##' @param seed if not NULL, this seed for the pseudo-random number generator will be used to sample as many seeds as the number of replicates, these new seeds being used to shuffle genotypes before partitioning per fold
##' @param nb.cores.rep number of cores to launch replicates in parallel (via \code{\link[parallel]{mclapply}}, on Unix-like computers)
##' @param nb.folds number of folds
##' @param remove.files remove files per fold (none/some/all); use \code{"some"} in real-life applications in order to keep estimates of SNP effects per fold, thereby allowing to perform genomic prediction afterwards by averaging them
##' @param nb.cores.fold number of cores to launch folds in parallel (via \code{\link[parallel]{mclapply}}, on Unix-like computers, or \code{\link[parallel]{parLapply}} on Windows); you can also use \code{\link[parallel]{detectCores}}; will be set to 1 if \code{nb.cores.rep} is bigger than 1
##' @param cl object returned by \code{\link[parallel]{makeCluster}}, necessary only if \code{nb.cores.fold > 1} and the computer runs Windows; if NULL in such cases, will be created silently
##' @param verbose verbosity level (0/1/2); there will be a progress bar only for \code{verbose=1}
##' @return data frame
##' @author Timothee Flutre
##' @seealso \code{\link{crossValWithGs3}}
##' @export
crossValRepWithGs3 <- function(genos,
                               dat,
                               config,
                               task.id="GS3",
                               binary.trait=FALSE,
                               ped.file="",
                               afs=NULL,
                               nb.reps=50,
                               seed=NULL,
                               nb.cores.rep=1,
                               nb.folds=10,
                               remove.files="some",
                               nb.cores.fold=1,
                               cl=NULL,
                               verbose=1){
  requireNamespace("parallel")
  stopifnot(is.numeric(nb.reps),
            nb.reps > 0)
  stopifnot(! is.null(ped.file))
  if(is.na(ped.file))
    ped.file <- ""

  if(! is.null(seed))
    set.seed(seed)
  seeds <- sample.int(n=10^5, size=nb.reps)

  if(nb.cores.rep > 1)
    nb.cores.fold <- 1

  if(is.null(afs)){
    if(verbose > 0){
      msg <- "estimate allele frequencies..."
      write(msg, stdout())
    }
    afs <- colMeans(genos, na.rm=TRUE) / 2
    if(verbose > 0)
      print(summary(afs))
  } else{
    stopifnot(is.numeric(afs),
              ! is.null(names(afs)),
              all(colnames(genos) %in% names(afs)))
    afs <- afs[colnames(genos)]
    stopifnot(all(! is.na(afs)),
              all(afs >= 0),
              all(afs <= 1))
  }

  out <- parallel::mclapply(seq(nb.reps), function(rep.id){
    if(verbose > 0){
      msg <- sprintf(fmt=paste0("rep %0", floor(log10(nb.reps))+1, "i",
                                "/", nb.reps),
                     rep.id)
      write(msg, stdout())
    }
    crossValWithGs3(genos=genos,
                    dat=dat,
                    config=config,
                    task.id=task.id,
                    ped.file=ped.file,
                    afs=afs,
                    rep.id=rep.id,
                    nb.folds=nb.folds,
                    seed=seeds[rep.id],
                    remove.files=remove.files,
                    nb.cores=nb.cores.fold,
                    verbose=verbose - 1)
  }, mc.cores=nb.cores.rep,
  mc.silent=ifelse(verbose == 0, TRUE, FALSE))

  out <- as.data.frame(cbind(rep=rep(1:nb.reps, each=nb.folds),
                             do.call(rbind, out)))
  attr(out, "seeds") <- seeds

  return(out)
}
