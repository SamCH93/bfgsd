## set R version (https://hub.docker.com/r/rocker/verse/tags)
FROM rocker/verse:4.5.2

## set up directories
RUN mkdir /home/rstudio/paper

## install R packages from CRAN the last day of the specified R version
RUN install2.r --error --skipinstalled --ncpus -1 \
    remotes knitr ggplot2 dplyr ggpubr mvtnorm xtable rpact scales tidyr ggrain && \
    ## TODO install from GitHub with fixed commit hash
    R -e "remotes::install_github(repo = 'SamCH93/bfpwr', subdir = 'package', ref = 'gsd', upgrade_dependencies = FALSE)"
