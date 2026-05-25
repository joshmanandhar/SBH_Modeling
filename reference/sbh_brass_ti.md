# Compute t(i) for calculated q(x)'s

This function takes q(x)'s calculated from the Brass (Trussell Variant)
method and Coale_demeny life tables to generate a reference year for a
specific q(x) based on the census year. for more information on the
formula, refer to the demo vignette.

## Usage

``` r
sbh_brass_ti(region_dataframe, i_brass, family, census_year)
```

## Arguments

- region_dataframe:

  a dataframe for a single region/area with SBH information for 5.yr age
  bands & both sexes from 15-19 up to 45-49, and a corresponding index
  column (i) for the age groups... columns should be named exactly:
  region, i, num_women, total_both, dead_both. all less region should be
  numeric. the dataframe should be 7 rows.

- i_brass:

  i corresponding to the q(x) of interest for the reference year. for
  example, i = 1 (ages 15-19) corresponds to x = 1, i = 2 (ages 15-19)
  corresponds to x = 2, etc. this should be a single integer ranging
  from 1-7 (inclusive).

- family:

  string representing family of Coale-Demeny life tables to be used
  (options: East, West, North, South - case sensitive)

- census_year:

  single integer representing the year of the census being used

## Value

the reference year for a given i corresponding to the q(x) of interest

## Examples

``` r
# Example SBH data
sbh_example <- data.frame(
  region = rep("ExampleRegion", 7),
  i = 1:7,
  num_women = c(100, 110, 120, 130, 120, 100, 90),
  total_both = c(20, 80, 150, 220, 260, 240, 180),
  dead_both = c(2, 5, 9, 14, 18, 20, 17)
)

# Compute reference time t(i) for i = 2
sbh_brass_ti(
  region_dataframe = sbh_example,
  i_brass = 2,
  family = "West",
  census_year = 2020
)
#> [1] 2016.99
```
