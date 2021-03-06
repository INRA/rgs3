% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gs3.R
\name{vcs2mcmc}
\alias{vcs2mcmc}
\title{Load GS3 results}
\usage{
vcs2mcmc(config, afs = NULL)
}
\arguments{
\item{config}{list containing the configuration for GS3}

\item{afs}{vector of allele frequencies (of all the SNPs which genotypes were provided, and only them); used to compute the variance of additive (and dominance) genotypic values from the variance of additive (and dominance) SNP effects (see Vitezica et al, 2013), and to compute narrow-sense heritability}
}
\value{
\code{\link[coda]{mcmc.list}}
}
\description{
Read the file containing the variance components' samples into a \code{\link[coda]{mcmc.list}} object.
}
\examples{
\dontrun{vcs <- vcs2mcmc(config$vcs.file)
summary(vcs)
coda::effectiveSize(vcs)
genos <- as.matrix(read.table(genos.file))
afs <- colMeans(genos) / 2
vcs <- vcs2mcmc(config$vcs.file, afs)
summary(vcs)}
}
\seealso{
\code{\link{execGs3}}
}
\author{
Timothee Flutre
}
