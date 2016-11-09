library(rgs3)
context("Gs3")

test_that("writeDataForGs3", {
  tmpd <- tempdir()

  phenos <- data.frame(geno=paste0("geno", 1:3),
                       trait.cont=c(0.24, 1.17, NA),
                       trait.bin=c(1, NA, 0))
  phenos.file <- paste0(tmpd, "/phenos.txt")
  inds <- setNames(object=1:3, nm=paste0("geno", 1:3))

  expected <- as.data.frame(matrix(c(1:3,
                                     c(0.24, 1.17, -9999),
                                     c(1, 0, 0)),
                                   nrow=3, ncol=3))

  writeDataForGs3(x=phenos,
                  file=phenos.file,
                  inds=inds,
                  col.id=1,
                  col.traits=c(2,3),
                  binary.traits=c(FALSE,TRUE))

  observed <- read.table(file=phenos.file)

  expect_equal(observed, expected)

  if(file.exists(phenos.file))
    file.remove(phenos.file)
})
