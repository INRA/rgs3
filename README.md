# Package "rgs3"

This directory contains the `rgs3` package for the [R](https://en.wikipedia.org/wiki/R_(programming_language)) software environment.
This package wraps the GS3 program for genomic selection.
See the [web page](http://snp.toulouse.inra.fr/~alegarra/) from Andrés Legarra for more information.

The development is funded by the Institut National de la Recherche Agronomique ([INRA](http://inra.fr/en)).
The copyright hence is owned by the INRA.
See the COPYING file for usage permissions.

The content of this directory is versioned using [git](https://en.wikipedia.org/wiki/Git_(software)), the central repository being hosted [here](https://github.com/timflutre/rgs3), on [GitHub](https://en.wikipedia.org/wiki/GitHub).
You can report [issues](https://github.com/timflutre/rgs3/issues) online, but remember to copy-paste the output of `sessionInfo()`.

When using this package, the latest version of the GS3 program should already be installed on your computer.
More precisely, the executable should be present in your [PATH](https://en.wikipedia.org/wiki/PATH_%28variable%29), under the name `gs3.exe` for Windows, and `gs3` for other operating systems (GNU/Linux, Mac OS).

* On Windows, you can save the executable in a new directory named `GS3`, for instance in `C:\Program Files`, and then add the path to this new directory to the environment variable `Path` (go to `Configuration parameters -> System -> Advanced`, or something similar).

* On other operating systems, you can save the executable in a new directory named `bin` in your home, and then add the path to this new directory to the environment variable `PATH` (use your `~/.bash_profile`).

To check if R properly detects the new directory, open a new R session, and call `Sys.getenv("PATH")`.
To check if the executable is found in your PATH, open a new R session, and call `system("gs3")` (or `system("gs3.exe")` for Windows).

As a lot of time and effort were spent in creating the GS3 program, please cite it when using it for data analysis:
```
Legarra, A., Ricard, A., Filangi, O. GS3, a  software  for  genome-wide genetic evaluations and validations. 2014.
```

For users, to install the `rgs3` package, the easiest is to install it directly from GitHub:
```
R> library(devtools)
R> install_github("timflutre/rgs3", build_vignettes=TRUE)
```

Note that creating the vignettes may take a couple of minutes.

Once the package is installed, you can start using it:
```
R> library(rgs3)
R> help(package="rgs3")
R> browseVignettes("rgs3")
```

You can cite this R package:
```
R> citation("rgs3")
```

See also `citation()` for citing R itself.

For developpers, when editing the content of this repo, increment the version of the package in the `DESCRIPTION` file, and execute the following commands:
```
$ Rscript -e 'library(devtools); devtools::document()'
$ R CMD build rgs3 # add '-no-build-vignettes' if necessary
$ R CMD check rgs3_<version>.tar.gz # add '--no-vignettes --no-build-vignettes --ignore-vignettes' if necessary
$ R CMD INSTALL rgs3_<version>.tar.gz
```
