# Compute t(i), the reference year, for a calculated q(x)

This function similarly uses the Brass (Trussell Variant) method and
Coale-Demeny life tables to generate a reference year for a specific
q(x) based on the census year. This is because a given q(x) being
computed won't necessarily correspond to the year the census was taken.
Once again, the Coale-Demeny tables can be found in the UN Manual X:
Indirect Techniques for Demographic Estimation but are also embedded in
the function directly. For more information on the formula, refer to the
SBH_Demo vignette.

## Usage

``` r
sbh_brass_ti(region_dataframe, i_brass, family, census_year)
```

## Arguments

- region_dataframe:

  A dataframe for a single area (i.e. district) with SBH information for
  both sexes (of children) for women in 5-year age bands ranging 15-19
  up to 45-49, and a corresponding index column (i) for the age groups.
  Columns should be named exactly: region, i, num_women, total_both
  (referring to Children Ever Born), dead_both (referring to Children
  Dead). All less region should be numeric. The dataframe should be 7
  rows.

- i_brass:

  'i' corresponding to the q(x) being calculated. For example, i = 1
  (ages 15-19) corresponds to x = 1 in q(x), i = 2 (ages 15-19)
  corresponds to x = 2, etc. This should be a single integer ranging
  from 1-7 (inclusive).

- family:

  A string field representing the family of Coale-Demeny life tables to
  be used (options: East, West, North, South - note this is case
  sensitive).

- census_year:

  A single integer representing the year of the census being analyzed.

## Value

A single value representing the reference year for a given i
corresponding to the q(x) of interest.

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
