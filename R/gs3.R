##' Write data
##'
##' Write the data file for GS3.
##' @param x data.frame containing individual identifiers (names of parameter "inds"), trait values (possibly several traits), covariables and cross-classified factors. Missing data coded as "NA" will be replaced by "0" if it is in a 'trait' column and the trait is binary, and by "-9999" otherwise.
##' @param file path to the text file to which x will be written
##' @param inds named vector with values from 1 to I and names the corresponding individual identifiers
##' @param col.id column index of the individual identifiers in x
##' @param col.traits column indices of traits in x
##' @param binary.traits logical vector indicating if traits are binary or not (should be of the same length as col.traits)
##' @return nothing
##' @author Timothee Flutre
##' @export
writeDataForGs3 <- function(x, file, inds, col.id=1, col.traits=2, binary.traits=FALSE){
  stopifnot(is.data.frame(x),
            ncol(x) >= 2, # at least individual identifiers and trait values
            is.character(file),
            is.vector(inds),
            is.numeric(inds),
            ! is.null(names(inds)),
            all(inds == 1:length(inds)),
            length(col.id) == 1,
            col.id >= 1,
            col.id <= ncol(x),
            is.factor(x[,col.id]),
            all(levels(x[,col.id]) %in% names(inds)),
            length(col.traits) >= 1,
            length(col.traits) <= ncol(x) - 1,
            all(col.traits >= 1),
            all(col.traits <= ncol(x)),
            all(! is.factor(x[,col.traits])),
            ! col.id %in% col.traits,
            all(is.logical(binary.traits)),
            length(binary.traits) == length(col.traits))
  if(any(binary.traits))
    stopifnot(all(x[,col.traits] %in% c(0,1,2,NA)))

  ## handle individual identifiers
  x[,col.id] <- as.character(x[,col.id])
  x[,col.id] <- setNames(object=inds[match(x[,col.id], names(inds))],
                         nm=NULL)

  ## handle missing values
  for(c in 1:ncol(x)){
    if(c == col.id)
      next
    idx <- which(is.na(x[,c]))
    if(length(idx) > 0){
      if(c %in% col.traits){
        if(binary.traits[c]){
          x[idx, c] <- 0
        } else
          x[idx, c] <- -9999
      } else
        x[idx, c] <- -9999
    }
  }

  ## handle factor levels
  for(c in 1:ncol(x))
    if(c != col.id & ! c %in% col.traits & is.factor(x[,c]))
      levels(x[,c]) <- 1:nlevels(x[,c])

  write.table(x=x,
              file=file,
              quote=FALSE,
              sep=" ",
              row.names=FALSE,
              col.names=FALSE)
}

##' Write genotypes
##'
##' Write the genotype file for GS3.
##' @param x matrix with SNP genotypes encoded as allele dose (i.e. 0/1/2 or NA/5), with individuals in rows and SNPs in columns. Row names are compulsory. Missing data coded as "NA" will be replaced by "5". Missing data already coded as "5" won't be modified.
##' @param file path to the text file to which x will be written
##' @param inds named vector with values from 1 to I and names the corresponding individual identifiers
##' @return nothing
##' @author Timothee Flutre
##' @export
writeGenosForGs3 <- function(x, file, inds){
  stopifnot(is.matrix(x),
            all(x %in% c(0,1,2,NA,5)),
            ! is.null(rownames(x)),
            is.vector(inds),
            is.numeric(inds),
            ! is.null(names(inds)),
            all(inds == 1:length(inds)),
            all(rownames(x) %in% names(inds)))

  ## handle the individual identifiers for GS3
  rownames(x) <- setNames(object=inds[match(rownames(x), names(inds))],
                          nm=NULL)

  ## handle missing genotypes for GS3
  x[is.na(x)] <- 5

  ## reformat (a bit slow)
  tmp <- cbind(sprintf("% 47s", rownames(x)),
               t(t(apply(x, 1, paste0, collapse=""))))

  write.table(x=tmp,
              file=file,
              quote=FALSE,
              sep=" ",
              row.names=FALSE,
              col.names=FALSE)
}

##' Write configuration
##'
##' Write the configuration file for GS3.
##' @param config.file path to the text file to which the configuration for GS3 will be written
##' @param data.file path to the text file with the data
##' @param ped.file path to the text file with the pedigree
##' @param genos.file path to the text file with the genotypes
##' @param num.loci number of loci
##' @param method BLUP/MCMCBLUP/VCE/PREDICT
##' @param simul one-letter character indicating if simulations should be performed (T) or not (F)
##' @param niter number of iterations for the Gibbs sampler
##' @param burnin burn-in for the Gibbs sampler
##' @param thin thinning for the Gibbs sampler
##' @param conv.crit convergence criterion (meaningful if BLUP)
##' @param correct correction (to avoid numerical problems)
##' @param vcs.file path to the text file to which the variance component samples will be written
##' @param sol.file path to the text file to which the solutions will be written
##' @param twc 2-element vector which first element corresponds to the column index of the trait values in the data file, and the second to the column index of the weights in the data file (use 0 if no weight)
##' @param num.eff number of effects
##' @param ptl 3-column data.frame indicating, for each covariable/factor, the position in the data file, type of effect and number of levels
##' @param vc 3-column data.frame indicating, for each variance component, the expected value and degrees of freedom
##' @param rec.id 1-element vector with a unique number for each record. Used  to  trace  the  records  across  the  cross-validation process.
##' @param cont vector with T or F indicating if the MCMC run is a continuation of a previous, interrupted one.
##' @param mod vector with T or F for each covariable/factor indicating if it has to be included or not in the model
##' @param ap prior proportions of the BayesCPi mixture
##' @param dp prior proportions of the BayesCPi mixture
##' @param use.mix one-letter character indicating if the Bayes C pi prior should be used (T) or not (F)
##' @param blasso boolean indicating if the Bayesian lasso prior should be used or not
##' @return nothing
##' @author Timothee Flutre
##' @export
writeConfigForGs3 <- function(config.file, data.file, ped.file=NULL, genos.file,
                              num.loci, method, simul="F",
                              niter=10000, burnin=2000, thin=10,
                              conv.crit="1d-8", correct=1000,
                              vcs.file="var.txt", sol.file="sol.txt",
                              twc, num.eff, ptl,
                              vc=data.frame(var=c("vara","vard","varg","varp","vare"),
                                            exp=c("2.52d-04","1.75d-06","3.56","2.15","0.19"),
                                            df=rep("-2", 5),
                                            stringsAsFactors=FALSE),
                              rec.id, cont="F", mod,
                              ap=c(1,10), dp=c(1,1), use.mix="F",
                              blasso=FALSE){
  stopifnot(file.exists(data.file),
            file.exists(genos.file))
  if(! is.null(ped.file))
    stopifnot(file.exists(ped.file))
  stopifnot(method %in% c("BLUP", "MCMCBLUP", "VCE", "PREDICT"),
            is.vector(twc),
            length(twc) == 2,
            is.data.frame(ptl),
            ncol(ptl) == 3,
            all(colnames(ptl) %in% c("position", "type", "nlevels")),
            all(ptl$type %in% c("cross", "cov", "add_SNP", "dom_SNP", "add_animal", "perm_diagonal")),
            nrow(ptl) == num.eff,
            is.vector(mod),
            length(mod) == num.eff,
            all(mod %in% c("T", "F")),
            is.data.frame(vc),
            ncol(vc) == 3,
            all(colnames(vc) %in% c("var", "exp", "df")),
            is.vector(ap),
            length(ap) == 2,
            is.vector(dp),
            length(dp) == 2,
            is.character(use.mix),
            length(use.mix) == 1,
            use.mix %in% c("T", "F"),
            is.logical(blasso),
            length(blasso) == 1)
  if(blasso)
    stopifnot(method == "VCE",
              use.mix == "F")

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
                "\n", num.loci)
  txt <- paste0(txt, "\nMETHOD (BLUP/MCMCBLUP/VCE/PREDICT)",
                "\n", method)
  txt <- paste0(txt, "\nSIMULATION",
                "\n", simul)
  txt <- paste0(txt, "\nGIBBS SAMPLING PARAMETERS")
  txt <- paste0(txt, "\nNITER",
                "\n", niter)
  txt <- paste0(txt, "\nBURNIN",
                "\n", burnin)
  txt <- paste0(txt, "\nTHIN",
                "\n", thin)
  txt <- paste0(txt, "\nCONV_CRIT (MEANINGFUL IF BLUP)",
                "\n", conv.crit)
  txt <- paste0(txt, "\nCORRECTION (to avoid numerical problems)",
                "\n", correct)
  txt <- paste0(txt, "\nVARIANCE COMPONENTS SAMPLES",
                "\n", vcs.file)
  txt <- paste0(txt, "\nSOLUTION FILE",
                "\n", sol.file)
  txt <- paste0(txt, "\nTRAIT AND WEIGHT COLUMNS",
                "\n", paste0(twc, collapse=" "))
  txt <- paste0(txt, "\nNUMBER OF EFFECTS",
                "\n", num.eff)
  txt <- paste0(txt, "\nPOSITION IN DATA FILE TYPE OF EFFECT  NUMBER OF LEVELS")
  for(i in 1:nrow(ptl))
    txt <- paste0(txt, "\n", ptl$position[i], " ", ptl$type[i], " ", ptl$nlevels[i])
  txt <- paste0(txt, "\nVARIANCE COMPONENTS (fixed for any BLUP, starting values for VCE)")
  for(i in 1:nrow(vc))
    txt <- paste0(txt, "\n", vc$var[i],
                  "\n", vc$exp[i], " ", vc$df[i])
  txt <- paste0(txt, "\nRECORD ID",
                "\n", rec.id)
  txt <- paste0(txt, "\nCONTINUATION (T/F)",
                "\n", cont)
  txt <- paste0(txt, "\nMODEL (T/F for each effect)",
                "\n", paste0(mod, collapse=" "))
  txt <- paste0(txt, "\nA PRIORI a",
                "\n", paste0(ap, collapse=" "))
  txt <- paste0(txt, "\na PRIORI D",
                "\n", paste0(dp, collapse=" "))
  txt <- paste0(txt, "\nUSE MIXTURE",
                "\n", use.mix)
  if(blasso)
    txt <- paste0(txt, "\n#OPTION BayesianLasso Tibshirani")
  txt <- paste0(txt, "\n")

  cat(txt, file=config.file)
}

##' Execute GS3
##'
##' Execute GS3 via a system call.
##' @param config.file path to the text file containing the configuration for GS3
##' @param stdouterr.file path to the text file to which the stdout and stderr will be written
##' @return return value (0 if success)
##' @author Timothee Flutre
##' @export
execGs3 <- function(config.file, stdouterr.file="gs3_stdouterr.txt"){
  stopifnot(file.exists(config.file))

  exec.name <- "gs3"
  inv <- TRUE
  if(Sys.info()["sysname"] == "Windows"){
    exec.name <- "gs3.exe"
    inv <- FALSE
  }

  ret <- system2(exec.name, args=c(config.file),
                 stdout=stdouterr.file, stderr=stdouterr.file, wait=TRUE,
                 invisible=inv)

  return(invisible(ret))
}

##' Load GS3 results
##'
##' Read the file containing the variance components' samples into a \code{\link[coda]{mcmc.list}} object.
##' @param vcs.file path to the file containing the variance components' samples
##' @return \code{\link[coda]{mcmc.list}}
##' @author Timothee Flutre
##' @examples
##' \dontrun{vcs <- vcs2mcmc(vcs.file)
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

  d <- read.table(vcs.file, header=TRUE)
  for(j in seq_along(d))
    if(! is.numeric(d[[j]]))
      d[[j]] <- NULL

  return(coda::mcmc.list(coda::mcmc(d)))
}
