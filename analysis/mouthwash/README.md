<!-- README.md is generated from README.Rmd. Please edit that file -->
Description
-----------

This Shiny app will let you run MOUTHWASH (Maximizing Over Unobservables To Help With Adaptive SHrinkage) interactively.

Running the Shiny App
---------------------

First, install all dependencies by typing the following code in `R`:

``` r
install.packages(c("shiny", "devtools", "SQUAREM", "cate", "ggplot2"))
devtools::install_github("stephens999/ashr", ref = "uni")
devtools::install_github("NKweiwang/flash")
devtools::install_github("dcgerard/succotashr")
source("https://bioconductor.org/biocLite.R")
biocLite(c("limma", "sva"))
```

Finally, to run the Shiny app, run in `R`:

``` r
shiny::runGitHub(repo = "succotash_sims", username = "dcgerard", subdir = "analysis/mouthwash")
```
