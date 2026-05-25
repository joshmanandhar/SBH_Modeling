# SBHModeling

# Purpose

This is a package created to work on obtaining mortality rate estimates
(and eventually uncertainty as well) from Summary Birth History data to
be used in small area estimation as covariates. This is under the
supervision of Prof. Jon Wakefield for Josh Manandhar’s work.

- The URL to the GitHub (i.e., the source code) is:
  <https://github.com/joshmanandhar/SBH_Modeling>
- The URL to the Pkgdown webpage is:
  (<https://joshmanandhar.github.io/SBH_Modeling>)

# How to install

This package is called `SBHModeling`. To install, run the following code
(in R):

``` r

library(devtools)
devtools::install_github("joshmanandhar/SBHModeling")
```

Upon completion, you can run the following code (in R):

``` r

library(SBHModeling)
```

Examples/walkthroughs can be found in the vignettes folder (files with
DEMO in the name)

# Dependencies

The package depends on the following packages: `dplyr`, `tidyr`,
`tidyverse`, `ggplot2`, `readxl`, `patchwork`, `stats`

# Session info

This package was developed in the following environment

``` r
devtools::session_info()
─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.5.2 (2025-10-31)
 os       macOS Sonoma 14.3
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/New_York
 date     2026-05-24
 rstudio  2026.01.0+392 Apple Blossom (desktop)
 pandoc   3.6.3 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
 quarto   1.8.25 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package     * version date (UTC) lib source
 askpass       1.2.1   2024-10-04 [1] CRAN (R 4.5.0)
 cachem        1.1.0   2024-05-16 [1] CRAN (R 4.5.0)
 cli           3.6.6   2026-04-09 [1] CRAN (R 4.5.2)
 credentials   2.0.3   2025-09-12 [1] CRAN (R 4.5.0)
 curl          7.0.0   2025-08-19 [1] CRAN (R 4.5.0)
 desc          1.4.3   2023-12-10 [1] CRAN (R 4.5.0)
 devtools      2.5.2   2026-04-30 [1] CRAN (R 4.5.2)
 digest        0.6.39  2025-11-19 [1] CRAN (R 4.5.2)
 ellipsis      0.3.3   2026-04-04 [1] CRAN (R 4.5.2)
 evaluate      1.0.5   2025-08-27 [1] CRAN (R 4.5.0)
 fastmap       1.2.0   2024-05-15 [1] CRAN (R 4.5.0)
 fs            2.1.0   2026-04-18 [1] CRAN (R 4.5.2)
 gert          2.3.1   2026-01-11 [1] CRAN (R 4.5.2)
 gh            1.5.0   2025-05-26 [1] CRAN (R 4.5.0)
 gitcreds      0.1.2   2022-09-08 [1] CRAN (R 4.5.0)
 glue          1.8.0   2024-09-30 [1] CRAN (R 4.5.0)
 htmltools     0.5.9   2025-12-04 [1] CRAN (R 4.5.2)
 httr2         1.2.2   2025-12-08 [1] CRAN (R 4.5.2)
 jsonlite      2.0.0   2025-03-27 [1] CRAN (R 4.5.0)
 knitr         1.51    2025-12-20 [1] CRAN (R 4.5.2)
 lifecycle     1.0.5   2026-01-08 [1] CRAN (R 4.5.2)
 magrittr      2.0.4   2025-09-12 [1] CRAN (R 4.5.0)
 memoise       2.0.1   2021-11-26 [1] CRAN (R 4.5.0)
 openssl       2.3.4   2025-09-30 [1] CRAN (R 4.5.0)
 otel          0.2.0   2025-08-29 [1] CRAN (R 4.5.0)
 pillar        1.11.1  2025-09-17 [1] CRAN (R 4.5.0)
 pkgbuild      1.4.8   2025-05-26 [1] CRAN (R 4.5.0)
 pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.5.0)
 pkgload       1.5.2   2026-04-22 [1] CRAN (R 4.5.2)
 purrr         1.2.1   2026-01-09 [1] CRAN (R 4.5.2)
 R6            2.6.1   2025-02-15 [1] CRAN (R 4.5.0)
 rappdirs      0.3.4   2026-01-17 [1] CRAN (R 4.5.2)
 rlang         1.2.0   2026-04-06 [1] CRAN (R 4.5.2)
 rmarkdown     2.31    2026-03-26 [1] CRAN (R 4.5.2)
 rprojroot     2.1.1   2025-08-26 [1] CRAN (R 4.5.0)
 rstudioapi    0.18.0  2026-01-16 [1] CRAN (R 4.5.2)
 sessioninfo   1.2.3   2025-02-05 [1] CRAN (R 4.5.0)
 sys           3.4.3   2024-10-04 [1] CRAN (R 4.5.0)
 tibble        3.3.1   2026-01-11 [1] CRAN (R 4.5.2)
 usethis       3.2.1   2025-09-06 [1] CRAN (R 4.5.0)
 vctrs         0.7.1   2026-01-23 [1] CRAN (R 4.5.2)
 whisker       0.4.1   2022-12-05 [1] CRAN (R 4.5.0)
 withr         3.0.2   2024-10-28 [1] CRAN (R 4.5.0)
 xfun          0.56    2026-01-18 [1] CRAN (R 4.5.2)
 yaml          2.3.12  2025-12-10 [1] CRAN (R 4.5.2)

 [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
```
