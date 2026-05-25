# Conversion of q(x) estimates to a q(y) of interest

This function takes q(x) estimates for a given x and uses Coale-Demeny
life tables to convert it to a different q(y) of interest for that given
reference year. for more information on the process, refer to the demo
vignette.

## Usage

``` r
sbh_interpolation_complete(sbh_result_table, family, y)
```

## Arguments

- sbh_result_table:

  a table with previously estimated q(x)s and their reference years for
  all age groups/i for one specific region. columns should be named
  region, q_i, i, x, and ref_year exactly. dataframe should be 7 rows
  long.

- family:

  the family for the life death tables to be used (includes North,
  South, West, East - case sensitive)

- y:

  the under-y mortality rate of interest that all previous q(x)s will be
  mapped to. this is NOT the same as i. for example, for u5mr you. input
  y = 5 and not y = 4. similarly, for under-10 you just input y = 10,
  not y = 5.

## Value

the sbh_result_table is returned with a new appended column representing
the interpolated mortality rates.

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

# Compute Brass q(x) estimates and reference years
region_name <- character(7)
ks <- numeric(7)
qs <- numeric(7)
yrs <- numeric(7)

for (i in 1:7) {
  result_brassq <- sbh_brass_qx(
    region_dataframe = sbh_example,
    i_brass = i,
    family = "West"
  )

  region_name[i] <- result_brassq$region
  ks[i] <- result_brassq$ki
  qs[i] <- result_brassq$qi

  yrs[i] <- sbh_brass_ti(
    region_dataframe = sbh_example,
    i_brass = i,
    family = "West",
    census_year = 2020
  )
}

# Create results table
sbh_results <- data.frame(
  region = region_name,
  q_i = qs,
  k_i = ks,
  i = 1:7,
  x = c(1, 2, 3, 5, 10, 15, 20),
  ref_year = yrs
)

# Interpolate all estimates to under-5 mortality
sbh_interpolation_complete(
  sbh_result_table = sbh_results,
  family = "West",
  y = 5
)
#> # A tibble: 7 × 7
#> # Rowwise: 
#>   region           q_i   k_i     i     x ref_year  q5_mr
#>   <chr>          <dbl> <dbl> <int> <dbl>    <dbl>  <dbl>
#> 1 ExampleRegion 0.0843 0.843     1     1    2019. 0.119 
#> 2 ExampleRegion 0.0597 0.955     2     2    2017. 0.0684
#> 3 ExampleRegion 0.0576 0.959     3     3    2015. 0.0618
#> 4 ExampleRegion 0.0629 0.988     4     5    2013. 0.0629
#> 5 ExampleRegion 0.0701 1.01      5    10    2010. 0.0642
#> 6 ExampleRegion 0.0835 1.00      6    15    2008. 0.0719
#> 7 ExampleRegion 0.0938 0.994     7    20    2005. 0.0737
```
