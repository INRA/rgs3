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

##' @import parallel

.onAttach <- function(libname, pkgname) {
  if(! requireNamespace("utils", quietly=TRUE))
    stop("Pkg utils needed for this function to work. Please install it.",
         call.=FALSE)
  msg <- paste0("package '", pkgname,
                "' (version ", utils::packageVersion(pkgname), ")",
                " is loaded",
                "\nCopyright 2016,2017 Institut National de la Recherche Agronomique (INRA)",
                "\nLicense GNU GPL version 3 or later")
  packageStartupMessage(msg)

  ## check that GS3 is installed
  exec.name <- "gs3"
  if(Sys.info()["sysname"] == "Windows")
    exec.name <- "gs3.exe"
  if(! file.exists(Sys.which(exec.name))){
    warning(paste0("can't find '", exec.name, "' in your PATH"))
  } else{
    ## check its version
    min.gs3.version <- "2.6.1"
    ret <- suppressWarnings(system2(command="gs3", args=c("dummyfile"),
                                    stdout=TRUE, stderr=FALSE))
    gs3.version <- trimws(ret[grep("        GS3        ", ret) + 1])
    if(utils::compareVersion(min.gs3.version, gs3.version) == 1)
      warning(paste0("installed GS3 is ", gs3.version,
                     " but should be >= ", min.gs3.version))
  }
}
