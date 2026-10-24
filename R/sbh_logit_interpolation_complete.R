#x-values for index..
x_values <- c(1, 2, 3, 4, 5, 10, 15, 20, 25, 30,
              35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
              85, 90, 95, 100, 105, 110, 115, 120, 125, 130)

logitlx_e60_north <- data.frame(
  x = x_values,
  combined_lx = c(
    -1.33002, -1.22732, -1.16637, -1.12258, -1.08999,
    -1.00905, -0.96635, -0.91382, -0.85013, -0.78726,
    -0.72475, -0.66077, -0.59046, -0.51475, -0.42188,
    -0.31373, -0.17374, 0.00771, 0.24439, 0.55421,
    0.97306, 1.55348, 2.30601, 3.45648, 4.78349,
    6.18826, 7.54854, 8.83708, 10.06797, 11.26005
  )
)

logitlx_e60_south <- data.frame(
  x = x_values,
  combined_lx = c(
    -1.13932, -1.02707, -0.98051, -0.95669, -0.94337,
    -0.91098, -0.88900, -0.85999, -0.82016, -0.78007,
    -0.73710, -0.69021, -0.63449, -0.56821, -0.48240,
    -0.37369, -0.22956, -0.04070, 0.21970, 0.58007,
    1.09073, 1.83304, 2.88179, 4.21903, 5.64094,
    6.99245, 8.26012, 9.46988, 10.64482, 11.79975
  )
)

logitlx_e60_east <- data.frame(
  x = x_values,
  combined_lx = c(
    -1.17947, -1.10549, -1.07580, -1.05684, -1.04341,
    -1.00433, -0.97971, -0.94156, -0.89016, -0.84220,
    -0.79105, -0.73425, -0.66826, -0.58661, -0.48112,
    -0.34745, -0.17916, 0.03487, 0.31235, 0.67727,
    1.17952, 1.88942, 2.87571, 4.15171, 5.56270,
    6.94208, 8.24338, 9.48044, 10.67480, 11.84277
  )
)

logitlx_e60_west <- data.frame(
  x = x_values,
  combined_lx = c(
    -1.29201, -1.19938, -1.15924, -1.13385, -1.11432,
    -1.06485, -1.02931, -0.97893, -0.91634, -0.85459,
    -0.79123, -0.72291, -0.64615, -0.55659, -0.44667,
    -0.31289, -0.14575, 0.06016, 0.32200, 0.66202,
    1.12034, 1.75553, 2.63151, 3.79176, 5.15205,
    6.55364, 7.89822, 9.17315, 10.39503, 11.58165
  )
)


sbh_interpolation_logit_row <- function(qxe_result, x, logit_lx_table, y){

  # #pull relevant logit lx for the given x and family
  standard_x_lifetable <- logit_lx_table$combined_lx[logit_lx_table$x == x]
  standard_y_lifetable <- logit_lx_table$combined_lx[logit_lx_table$x == y]

  # avoid 0 and 1 of qxe_result as this will mess up the logit function
  eps <- 1e-10
  qxe_result <- pmin(pmax(qxe_result, eps), 1 - eps)

  #convert the qxe_result to 1/2(logit) as outlined in the iussp method
  qxe_logit_trans <- 0.5 * log(qxe_result / (1- qxe_result))
  #take difference of this transofrmation with the pulled standard
  q_alpha <- qxe_logit_trans  - standard_x_lifetable

  qy_est <- exp(2*(q_alpha + standard_y_lifetable)) / (1 + exp(2*(q_alpha + standard_y_lifetable)))
  #qy_est <- plogis(2 * (q_alpha + standard_y_lifetable))
  return (qy_est)

  }

#' Interpolation of q(x) estimates to a specific q(y) of interest using the IUSSP method
#'
#' This function takes estimated q(x) values for a given region and uses interpolation with the help of model life tables to convert them all to a specific q(y) of interest (for their respective reference years).
#' For example, if the q(x) estimates are obtained for a region and one is interested in reconstructing q(5) - U5MR over the reference years.
#' Specifically, this function uses the estimated probabilities' conversion to the value a: a level parameter together with logit transformations of standard life tables.
#' The logit-trnasformed life tables are those using a life expectancy of 60 years, as highlighted and linked in the IUSSP method, but are also embedded in the function for ease.
#' For more information on the process and formulas, refer to the SBH_Demo vignette.
#'
#' @param sbh_result_table A dataframe with previously constructed q(x) estimates and their reference years for all age groups (i's) for one specific region.
#' This can be obtained by using the sbh_brass_qx and sbh_brass_ti functions. Refer to the example below or the SBH_Demo vignette to view the process for obtaining such from SBH information.
#' Columns should be named region, q_i, i, x, and ref_year exactly, with all except region being numeric. The dataframe should be 7 rows long.
#'
#' @param family  A string field representing the family of Coale-Demeny life tables to be used (options: East, West, North, South - note this is case sensitive).
#' @param y The under-y mortality rate of interest that all previous q(x)s will be interpolated to.
#' This is NOT the same as i. For example, for U5MR (Under-5 Mortality) one should input y = 5 and not y = 4. Similarly, for under-10 mortality, you just input y = 10 and not y = 5.
#' Recall previously in the q(x) functions however, you would input the i that corresponds to such mortality age.
#'
#' @returns The sbh_result_table itself is returned with a new appended column representing the interpolated mortality rates for the y of interest.
#'
#' @examples
#' # Example SBH data
#' sbh_example <- data.frame(
#'   region = rep("ExampleRegion", 7),
#'   i = 1:7,
#'   num_women = c(100, 110, 120, 130, 120, 100, 90),
#'   total_both = c(20, 80, 150, 220, 260, 240, 180),
#'   dead_both = c(2, 5, 9, 14, 18, 20, 17)
#' )
#'
#' # Compute Brass q(x) estimates and reference years
#' region_name <- character(7)
#' ks <- numeric(7)
#' qs <- numeric(7)
#' yrs <- numeric(7)
#'
#' for (i in 1:7) {
#'   result_brassq <- sbh_brass_qx(
#'     region_dataframe = sbh_example,
#'     i_brass = i,
#'     family = "West"
#'   )
#'
#'   region_name[i] <- result_brassq$region
#'   ks[i] <- result_brassq$ki
#'   qs[i] <- result_brassq$qi
#'
#'   yrs[i] <- sbh_brass_ti(
#'     region_dataframe = sbh_example,
#'     i_brass = i,
#'     family = "West",
#'     census_year = 2020
#'   )
#' }
#'
#' # Create results table
#' sbh_results <- data.frame(
#'   region = region_name,
#'   q_i = qs,
#'   k_i = ks,
#'   i = 1:7,
#'   x = c(1, 2, 3, 5, 10, 15, 20),
#'   ref_year = yrs
#' )
#'
#' # Interpolate all estimates to under-5 mortality
#' sbh_logit_interpolation_complete(
#'   sbh_result_table = sbh_results,
#'   family = "West",
#'   y = 5
#' )
#'
#' @export
sbh_logit_interpolation_complete <- function(sbh_result_table, family, y) {
#check inputs are as expected...
stopifnot(
  is.data.frame(sbh_result_table),
  all(c("region", "q_i", "i", "x", "ref_year") %in%
        names(sbh_result_table)),
  nrow(sbh_result_table) == 7,

  is.numeric(sbh_result_table$i),
  all(sbh_result_table$i == 1:7),

  is.numeric(sbh_result_table$x),
  all(sbh_result_table$x %in% c(1,2,3,5,10,15,20)),

  is.numeric(sbh_result_table$q_i),
  is.numeric(sbh_result_table$ref_year),
  length(unique(sbh_result_table$region)) == 1,

  length(family) == 1,
  family %in% c("West", "North", "East", "South"),

  length(y) == 1,
  is.numeric(y),
  y %% 1 == 0,
  y %in% c(1,2,3,5,10,15,20)
)

#definte new column name and pull the correct life table based on inputs
new_col <- paste0("q", y, "_mr")
if (family == "West") {
  logit_lx_table <- logitlx_e60_west
} else if (family == "North") {
  logit_lx_table <- logitlx_e60_north
} else if (family == "South") {
  logit_lx_table <- logitlx_e60_south
} else if (family == "East") {
  logit_lx_table <- logitlx_e60_east
}

#for each q(x) in the summary table, compute the interpolated q(y)
sbh_result_table |> dplyr::rowwise() |> dplyr::mutate(
  !!new_col := sbh_interpolation_logit_row(
    qxe_result = q_i,
    x  = x,
    logit_lx_table = logit_lx_table,
    y = y
  )
)
}

