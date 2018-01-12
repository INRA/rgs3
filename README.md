# Package "rgs3"

This directory contains the `rgs3` package for the [R](https://en.wikipedia.org/wiki/R_(programming_language)) software environment.
This package wraps the GS3 program for genomic prediction and selection.
See the [web page](http://snp.toulouse.inra.fr/~alegarra) from Andrés Legarra for more information.

The `rgs3` package is available under a free software license, the [GNU Public License](https://www.gnu.org/licenses/gpl.html) (version 3 or later).
See the COPYING file for details.
The copyright is owned by the [INRA](http://www.inra.fr).

The content of this directory is versioned using [git](https://en.wikipedia.org/wiki/Git_(software)), the central repository being hosted [here](https://github.com/INRA/rgs3), on [GitHub](https://en.wikipedia.org/wiki/GitHub), and the institutional repository being hosted [there](https://sourcesup.renater.fr/projects/rgs3/), on [SourceSup](https://sourcesup.renater.fr/).


## Installation

Before installing the `rgs3` package, the latest version of the GS3 program should already be installed on your computer (see on [GitHub](https://github.com/alegarra/gs3/releases/latest)).
More precisely, the executable should be present in your [PATH](https://en.wikipedia.org/wiki/PATH_%28variable%29), under the name `gs3` for Unix-like operating systems (GNU/Linux, Mac OS) and `gs3.exe` for Microsoft Windows.

* On Unix-like operating systems, you can save the executable in a new directory named `bin` in your [home directory](https://en.wikipedia.org/wiki/Home_directory), and then add the path to this new directory to the environment variable `PATH` (use your `~/.bash_profile`).

* On Windows, you can save the executable in a new directory named `GS3`, for instance in `C:\Program Files`, and then add the path to this new directory to the environment variable `Path` (go to `Configuration parameters -> System -> Advanced`, or something similar).

To check if R properly detects the new directory, open a new R session, and call `Sys.getenv("PATH")`.
To check if the executable is found in your PATH, open a new R session, and call `system("gs3")` (or `system("gs3.exe")` for Windows).

Then, to install the `rgs3` package, the easiest is to install it directly from GitHub.
Open an R session and run the following commands:
```
library(devtools) # can be installed from the CRAN
install_github("INRA/rgs3", build_vignettes=TRUE)
```

Note that creating the vignettes may take a couple of minutes.

Once this is done, the `rgs3` package should be available on your computer.


## Usage

Once the `rgs3` package is installed on your computer, it can be loaded into a R session:
```
library(rgs3)
help(package="rgs3")
browseVignettes("rgs3")
```


## Citation

As a lot of time and effort were spent in creating the GS3 program, please cite it when using it for data analysis:
```
Legarra, A., Ricard, A., Filangi, O. GS3, a  software  for  genome-wide genetic evaluations and validations. 2014.
```

You should also cite the `rgs3` package:
```
citation("rgs3")
```

See also `citation()` for citing R itself.


## Issues

When encountering a problem with the package, you can report issues on GitHub directly ([here](https://github.com/INRA/rgs3/issues)).
Remember to copy-paste the output of ` sessionInfo()` to help efficiently diagnose the problem and find a solution.


## Contributing

You can contribute in various ways:

* report an issue (online, see the above section);

* suggest improvements (in the same way as issues);

* propose a [pull request](https://github.com/INRA/rgs3/pulls) (after creating a new [branch](https://www.git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell)).
