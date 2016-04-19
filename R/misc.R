.onAttach <- function(libname, pkgname) {
  if(! requireNamespace("utils", quietly=TRUE))
    stop("Pkg utils needed for this function to work. Please install it.",
         call.=FALSE)
  msg <- paste0("package '", pkgname,
                "' (version ", utils::packageVersion(pkgname), ")",
                " is loaded",
                "\ndev at https://github.com/timflutre/rgs3")
  packageStartupMessage(msg)

  exec.name <- "gs3"
  if(Sys.info()["sysname"] == "Windows")
    exec.name <- "gs3.exe"
  if(! file.exists(Sys.which(exec.name)))
    warning(paste0("can't find '", exec.name, "' in your PATH"))
}
