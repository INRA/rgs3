## Copyright 2016,2017 Institut National de la Recherche Agronomique (INRA)
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
##' @param x data.frame containing individual identifiers (names of parameter "inds"), trait values (possibly several traits), covariables and cross-classified factors. Missing data coded as "NA" will be replaced by "0" if it is in a 'trait' column and the trait is binary, and by "-9999" otherwise.
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
getDefaultConfig <- function(nb.snps=NA, method="VCE", ptl=NULL, twc=NA,
                             rec.id=NA, use.mix="F", blasso=FALSE,
                             task.id="GS3"){
  if(! is.na(nb.snps))
    stopifnot(is.numeric(nb.snps),
              length(nb.snps) == 1,
              nb.snps > 0)
  if(! is.null(ptl)){
    stopifnot(is.data.frame(ptl),
              ncol(ptl) == 3,
              all(colnames(ptl) %in% c("position", "type", "nlevels")),
              all(ptl$type %in% c("cross", "cov", "add_SNP", "dom_SNP",
                                  "add_animal", "perm_diagonal")))
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

  config <- list(num.loci=nb.snps,
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
    if(all(c("num.loci", "method", "simul", "niter", "burnin", "thin",
             "conv.crit", "correct", "vcs.file", "sol.file", "twc", "num.eff",
             "ptl", "vc", "rec.id", "cont", "mod", "ap", "dp", "use.mix",
             "blasso") %in% names(config))){
      all(! is.na(config$num.loci),
          config$method %in% c("BLUP", "MCMCBLUP", "VCE", "PREDICT"),
          is.vector(config$twc),
          length(config$twc) == 2,
          config$twc[1] > 0,
          config$twc[2] >= 0,
          config$twc[1] != config$twc[2],
          is.data.frame(config$ptl),
          ncol(config$ptl) == 3,
          all(colnames(config$ptl) %in% c("position", "type", "nlevels")),
          all(config$ptl$type %in% c("cross", "cov", "add_SNP", "dom_SNP",
                                     "add_animal", "perm_diagonal")),
          nrow(config$ptl) == config$num.eff,
          is.vector(config$mod),
          length(config$mod) == config$num.eff,
          all(config$mod %in% c("T", "F")),
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
##'   \item{ptl}{3-column data.frame indicating, for each covariable/factor, its position (column) in the data file, type of effect and number of levels}
##'   \item{vc}{3-column data.frame indicating, for each variance component, its expected value and degrees of freedom}
##'   \item{rec.id}{1-element vector indicating the column in the data file corresponding to the genotype identifiers}
##'   \item{cont}{vector with T or F indicating if the MCMC run is a continuation of a previous, interrupted one}
##'   \item{mod}{vector with T or F for each covariable/factor indicating if it has to be included or not in the model}
##'   \item{ap}{prior proportions of the BayesCPi mixture}
##'   \item{dp}{prior proportions of the BayesCPi mixture}
##'   \item{use.mix}{one-letter character indicating if the Bayes C pi prior should be used (T) or not (F)}
##'   \item{blasso}{boolean indicating if the Bayesian lasso prior should be used or not}
##' }
##' @param data.file path to the text file with the data
##' @param ped.file path to the text file with the pedigree
##' @param genos.file path to the text file with the genotypes
##' @param task.id character containing the task identifier used as prefix for the configuration file; will be followed by \code{"_config.txt"}
##' @return path to the text file to which the configuration for GS3 will be written
##' @author Timothee Flutre
##' @seealso \code{\link{getDefaultConfig}}, \code{\link{execGs3}}
##' @export
writeConfigForGs3 <- function(config,
                              data.file,
                              ped.file=NULL,
                              genos.file=NULL,
                              task.id="GS3"){
  stopifnot(isValidConfig(config),
            file.exists(data.file))
  tmp <- utils::read.table(file=data.file, sep=" ", nrows=2)
  stopifnot(config$twc[1] <= ncol(tmp),
            ifelse(config$twc[2] != 0,
                   config$twc[2] <= ncol(tmp),
                   TRUE))
  if(! is.null(ped.file))
    stopifnot(file.exists(ped.file))
  if(! is.null(genos.file)){
    stopifnot(file.exists(genos.file))
    tmp <- readLines(con=genos.file, n=1)
    stopifnot(config$num.loci == nchar(tmp) - 48) # see writeGenosForGs3()
  }

  config.file <- paste0(task.id, "_config.txt")

  txt <- paste0("DATAFILE",
                "\n", data.file)
  if(! is.null(ped.file)){
    txt <- paste0(txt, "\nPEDIGREE FILE",
                  "\n", ped.file)
  } else
    txt <- paste0(txt, "\nPEDIGREE FILE",
                  "\n", "")
  txt <- paste0(txt, "\nGENOTYPE FILE",
                "\n", genos.file)
  txt <- paste0(txt, "\nNUMBER OF LOCI (might be 0)",
                "\n", config$num.loci)
  txt <- paste0(txt, "\nMETHOD (BLUP/MCMCBLUP/VCE/PREDICT)",
                "\n", config$method)
  txt <- paste0(txt, "\nSIMULATION",
                "\n", config$simul)
  txt <- paste0(txt, "\nGIBBS SAMPLING PARAMETERS")
  txt <- paste0(txt, "\nNITER",
                "\n", config$niter)
  txt <- paste0(txt, "\nBURNIN",
                "\n", config$burnin)
  txt <- paste0(txt, "\nTHIN",
                "\n", config$thin)
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

  if(file.exists("freq")){
    file.rename(from="freq",
                to=paste0(task.id, "_freq.txt"))
  }

  return(stdouterr.file)
}

##' Load GS3 results
##'
##' Read the file containing the variance components' samples into a \code{\link[coda]{mcmc.list}} object.
##' @param vcs.file path to the file containing the variance components' samples
##' @return \code{\link[coda]{mcmc.list}}
##' @author Timothee Flutre
##' @seealso \code{link{execGs3}}
##' @examples
##' \dontrun{vcs <- vcs2mcmc(config$vcs.file)
##' summary(vcs)
##' coda::effectiveSize(vcs)
##' genos <- as.matrix(read.table(genos.file))
##' afs <- colMeans(genos) / 2
##' tmp <- cbind(as.matrix(vcs[[1]]), varA=vcs[[1]][,"vara"] * 2 * sum(afs * (1 - afs)))
##' vcs[[1]] <- coda::as.mcmc(tmp)
##' summary(vcs)}
##' @export
vcs2mcmc <- function(vcs.file){
  if(! requireNamespace("coda", quietly=TRUE))
    stop("Pkg coda needed for this function to work. Please install it.",
         call.=FALSE)
  stopifnot(is.character(vcs.file),
            length(vcs.file) == 1,
            file.exists(vcs.file))

  d <- utils::read.table(vcs.file, header=TRUE, check.names=FALSE)
  for(j in seq_along(d))
    if(! is.numeric(d[[j]]))
      d[[j]] <- NULL

  return(coda::mcmc.list(coda::mcmc(d)))
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
                   paste0(config.file, "_EBVs"),
                   paste0(config.file, "_cont"),
                   paste0(config.file, "_finalEstimates"))
  files.tormv <- c(files.tormv, c(paste0(task.id, "_freq.txt"),
                                  paste0(task.id, "_stdouterr.txt")))

  ## remove the files
  for(f in files.tormv){
    if(file.exists(f))
      file.remove(f)
    else{
      msg <- paste0(" file '", f, "' doesn't exist")
      warning(msg)
    }
  }
}

##' Partition for cross-validation
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
##' @param remove.files remove files per fold (none/some/all); use \code{"some"} in real-life applications in order to keep estimates of SNP effects per fold, thereby allowing to perform genomic prediction afterwards by averaging them
##' @param verbose verbosity level (0/1/2); there will be a progress bar only for \code{verbose=1}
##' @return named vector with metrics
##' @author Timothee Flutre
##' @seealso \code{\link{crossValWithGs3}}
##' @export
crossValFold <- function(task.id, fold.ids, fold.id,
                         valid.geno.idx.per.fold,
                         dat, col.id, col.trait, binary.trait,
                         genos, config, ped.file,
                         remove.files, verbose){
  out <- stats::setNames(c(fold.id, rep(NA, 7)),
                         c("fold", "train.size", "valid.size",
                           "rmspe", "cor.p", "cor.s",
                           "reg.intercept", "reg.slope"))

  nb.folds <- length(fold.ids)

  ## set up file names
  tif <- paste0(task.id, "-", fold.ids[fold.id])
  dat.train.file <- paste0(tif, "_dat-train.txt")
  geno.train.file <- paste0(tif, "_geno-train.txt")
  dat.valid.file <- paste0(tif, "_dat-valid.txt")
  geno.valid.file <- paste0(tif, "_geno-valid.txt")

  ## training
  if(verbose > 1)
    write(paste0("fold ", fold.id, "/", nb.folds,
                 ": start training"), stdout())
  train.genos <- names(which(valid.geno.idx.per.fold != fold.id))
  out["train.size"] <- length(train.genos)
  dat.train <- droplevels(dat[dat[,col.id] %in% train.genos,])
  genos.train <- genos[train.genos,]
  inds.train <- stats::setNames(object=1:nlevels(dat.train[, col.id]),
                                nm=levels(dat.train[, col.id]))
  if(verbose > 1)
    write(paste0("fold ", fold.id, "/", nb.folds, ":",
                 " nb.genos=", nrow(genos.train),
                 " nb.snps=", ncol(genos.train),
                 " nb.obs=", nrow(dat.train)),
          stdout())
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
  config.train$vcs.file <- sub(task.id, tif, config$vcs.file)
  config.train$sol.file <- sub(task.id, tif, config$sol.file)
  config.train.file <- writeConfigForGs3(config = config.train,
                                         data.file = dat.train.file,
                                         ped.file = ped.file,
                                         genos.file = geno.train.file,
                                         task.id = tif)
  stdouterr.train.file <- execGs3(config.file = config.train.file,
                                  task.id = tif)
  if(FALSE) # for debugging purposes
    stdouterr.train <- readLines(stdouterr.train.file)
  if(verbose > 1)
    write(paste0("fold ", fold.id, "/", nb.folds,
                 ": end training"), stdout())

  ## validation
  if(verbose > 1)
    write(paste0("fold ", fold.id, "/", nb.folds,
                 ": start validation"), stdout())
  valid.genos <- names(which(valid.geno.idx.per.fold == fold.id))
  stopifnot(all(! valid.genos %in% train.genos))
  out["valid.size"] <- length(valid.genos)
  dat.valid <- droplevels(dat[dat[,col.id] %in% valid.genos,])
  genos.valid <- genos[valid.genos,]
  inds.valid <- stats::setNames(object=1:nlevels(dat.valid[, col.id]),
                                nm=levels(dat.valid[, col.id]))
  dat.valid.to.predict <- dat.valid
  dat.valid.to.predict[, col.trait] <- NA
  if(verbose > 1)
    write(paste0("fold ", fold.id, "/", nb.folds, ":",
                 " nb.genos=", nrow(genos.valid),
                 " nb.snps=", ncol(genos.valid),
                 " nb.obs=", nrow(dat.valid)),
          stdout())
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
  config.valid$method <- "PREDICT"
  config.valid.file <- writeConfigForGs3(config = config.valid,
                                         data.file = dat.valid.file,
                                         ped.file = ped.file,
                                         genos.file = geno.valid.file,
                                         task.id = tif)
  stdouterr.valid.file <- execGs3(config.file = config.valid.file,
                                  task.id = tif)
  if(FALSE) # for debugging purposes
    stdouterr.valid <- readLines(stdouterr.valid.file)
  pred.file <- paste0(tif, "_pred.txt")
  file.rename(from="predictions", to=pred.file)
  if(verbose > 1)
    write(paste0("fold ", fold.id, "/", nb.folds,
                 ": end validation"), stdout())

  ## assess accuracy
  if(FALSE){ # for debugging purposes
    vcs <- utils::read.table(file=config.train$vcs.file, header=TRUE,
                             check.names=FALSE)
    ebvs <- utils::read.table(file=paste0(config.valid.file, "_EBVs"),
                              header=TRUE)
    sols <- utils::read.table(file=config.train$sol.file, header=TRUE)
  }
  pred <- utils::read.table(file=pred.file, header=TRUE)
  stopifnot(nrow(pred) == nrow(dat.valid))
  pred$geno <- names(inds.valid)[match(pred$id, inds.valid)]
  pred <- pred[order(pred$geno),]
  dat.valid <- dat.valid[order(dat.valid[,col.id]),]
  stopifnot(all(pred$geno == as.character(dat.valid[,col.id])))
  out["rmspe"] <- sqrt(mean((dat.valid[, col.trait] -
                             pred$prediction)^2))
  out["cor.p"] <- stats::cor(x=dat.valid[, col.trait],
                             y=pred$prediction,
                             method="pearson")
  out["cor.s"] <- stats::cor(x=dat.valid[, col.trait],
                             y=pred$prediction,
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
                       paste0(config.train.file, "_EBVs"),
                       paste0(config.train.file, "_cont"),
                       dat.valid.file,
                       geno.valid.file,
                       config.valid.file,
                       paste0(config.valid.file, "_EBVs"),
                       stdouterr.valid.file,
                       pred.file,
                       paste0(tif, "_freq.txt"))
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

##' Cross-validation
##'
##' Perform K-fold cross-validation with GS3, and report metrics as advised in \href{http://dx.doi.org/10.1534/genetics.112.147983}{Daetwyler et al. (2013)}.
##' Files are saved in the current directory.
##' @param genos matrix of SNP genotypes
##' @param dat data frame with phenotypes
##' @param config list containing the configuration for GS3
##' @param task.id character containing the task identifier used as prefix for the output files (for each fold, its index will be added)
##' @param binary.trait logical
##' @param ped.file path to the file containing the pedigree
##' @param nb.folds number of folds
##' @param seed if not NULL, this seed for the pseudo-random number generator will be used to shuffle genotypes before partitioning per fold
##' @param remove.files remove files per fold (none/some/all); use \code{"some"} in real-life applications in order to keep estimates of SNP effects per fold, thereby allowing to perform genomic prediction afterwards by averaging them
##' @param nb.cores number of cores to launch folds in parallel (via \code{\link[parallel]{mclapply}}, on Unix-like computers, or \code{\link[parallel]{parLapply}} on Windows); you can also use \code{\link[parallel]{detectCores}}
##' @param cl object returned by \code{\link[parallel]{makeCluster}}, necessary only if \code{nb.cores > 1} and the computer runs Windows; if NULL in such cases, will be created silently
##' @param verbose verbosity level (0/1/2); there will be a progress bar only for \code{verbose=1}
##' @return data.frame
##' @author Helene Muranty, Timothee Flutre
##' @seealso \code{\link{getDefaultConfig}}, \code{\link{writeConfigForGs3}}, \code{\link{execGs3}}, \code{\link{getPartitionGenos}}
##' @export
crossValWithGs3 <- function(genos,
                            dat,
                            config,
                            task.id="GS3",
                            binary.trait=FALSE,
                            ped.file=NULL,
                            nb.folds=10,
                            seed=NULL,
                            remove.files="some",
                            nb.cores=1,
                            cl=NULL,
                            verbose=1){
  requireNamespace("parallel")
  stopifnot(isValidConfig(config))
  col.id <- config$rec.id
  col.trait <- config$twc[1]
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

  ## prepare output
  out <- data.frame(fold=1:nb.folds,
                    train.size=NA,
                    valid.size=NA,
                    ## rel=NA,
                    ## rel.sq=NA,
                    ## rel.top10=NA,
                    rmspe=NA,
                    cor.p=NA,
                    cor.s=NA,
                    reg.intercept=NA,
                    reg.slope=NA,
                    stringsAsFactors=FALSE)

  ## prepare partitions
  valid.geno.idx.per.fold <- getPartitionGenos(geno.names=rownames(genos),
                                               nb.folds=nb.folds,
                                               seed=seed)
  stopifnot(! any(duplicated(names(valid.geno.idx.per.fold))))

  ## perform cross-validation
  fold.ids <- sprintf(fmt=paste0("%0", floor(log10(nb.folds))+1, "i"),
                      1:nb.folds)
  ## out <- lapply(seq(nb.folds), function(fold.id){
  ##   crossValFold(task.id, fold.ids, fold.id,
  ##                valid.geno.idx.per.fold,
  ##                dat, col.id, col.trait, binary.trait,
  ##                genos, config, ped.file,
  ##                remove.files, verbose)
  ## })
  if(Sys.info()["sysname"] == "Windows"){
    out <- parallel::parLapply(cl, seq(nb.folds), function(fold.id){
      crossValFold(task.id, fold.ids, fold.id,
                   valid.geno.idx.per.fold,
                   dat, col.id, col.trait, binary.trait,
                   genos, config, ped.file,
                   remove.files, verbose)
    })
  } else{
    out <- parallel::mclapply(seq(nb.folds), function(fold.id){
      crossValFold(task.id, fold.ids, fold.id,
                   valid.geno.idx.per.fold,
                   dat, col.id, col.trait, binary.trait,
                   genos, config, ped.file,
                   remove.files, verbose)
    }, mc.cores=nb.cores)
  }
  out <- do.call(rbind, out)

  if(wasClCreated)
    parallel::stopCluster(cl)

  return(out)
}
