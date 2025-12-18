# Bayes Factor Group Sequential Designs

This repository contains code and data to reproduce result from the manuscript

> Pawel, S., Held. L. (2025). Bayes Factor Group Sequential Designs. <https://github.com/SamCH93/bfgsd>

To cite our work, use the following BibTeX reference

```BibTeX
@article{PawelHeld2025,
  year = {2025},
  author = {Samuel Pawel and Leonhard Held},
  title = {Bayes Factor Group Sequential Designs},
  url = {Bayes Factor Group Sequential Designs}
}
```

## Reproducing the results with Docker

Make sure to have Docker and Make installed, then run `make docker-rstudio` from
the root directory of this git repository. This will install all necessary
dependencies. RStudio Server can then be opened from a browser
(http://localhost:8787), and the R scripts in `/paper`, for example,
`/paper/BFGSD.R`, which contains all code for the results from the paper), can
be rerun. Make sure to change the working directory to `/paper` inside RStudio
Server before running the R scripts. Running `make docker-paper` produces the
`paper/BFGSD.tex` file from the `paper/BFGSD.Rnw` source file (dynamically
inserting numbers and figures) and then compiles it to a PDF (requires a local
LaTeX installation; only tested with TeX Live 2023/Debian).

## Reproducing the results locally

Make sure to have [R](https://www.r-project.org/) installed and then install the
required R packages with the R commands

``` r
install.packages(c("remotes", "knitr", "ggplot2", "dplyr", "ggpubr", "mvtnorm",
                   "xtable", "rpact", "scales", "tidyr", "ggrain"))
remotes::install_github(repo = "SamCH93/bfpwr", subdir = "package", ref = "gsd")
```

Go to the `/paper` diretory an run either the R script `/paper/BFGSD.R`, which
contains all code for the results from the paper, or knit and compile the
manuscript directly to a PDF (requires a local LaTeX installation) with the R
command

``` r
knitr::knit2pdf("BFGSD.Rnw")
```

This should reproduce the file `paper/BFGSD.pdf`
