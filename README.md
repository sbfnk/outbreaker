[![Travis-CI Build Status](https://travis-ci.org/thibautjombart/outbreaker.svg?branch=master)](https://travis-ci.org/thibautjombart/outbreaker)

# outbreaker
*Disease outbreak reconstruction from epidemiological and genetic data*

Welcome to the project page of *outbreaker*!
This page is dedicated to the development of the R package *outbreaker*.
Whenever possible, it is recommended to use the development version of the package rather than the older CRAN version. The devel version (at least the master branch) is meant to be as functional as the CRAN version, but integrates the latest new features and bug fixes.

For documentation and tutorials, see [The R-epi-project](https://sites.google.com/site/therepiproject/r-pac/outbreaker).


Installing *outbreaker* devel
-------------
To install the development version from github (the package *devtools is required*):

```r
library(devtools)
install_github("thibautjombart/outbreaker")
```
Note that on Windows, a toolkit ([Rtools](https://cran.r-project.org/bin/windows/Rtools/)) needs to be installed separately for *devtools* to work. On every system, [*GNU scientific library*](http://www.gnu.org/software/gsl/) also needs to be installed for *outbreaker* to be compiled. 

Once installed, the package can be loaded using:

```r
library("outbreaker")
```

Asking a question
------------------
- general questions can be asked on the [R-epi forum](http://groups.google.com/forum/#!forum/r-epi)
- for bug reports, feature requests, contributions, use github's [issue system](https://github.com/thibautjombart/outbreaker/issues)

