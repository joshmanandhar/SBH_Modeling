lifedeath_levels_west <- data.frame(
  level = 1:25,
  e0 = c(
    19.0, 21.4, 23.9, 26.4, 28.8,
    31.3, 33.7, 36.2, 38.6, 41.1,
    43.5, 46.0, 48.5, 51.0, 53.4,
    55.8, 58.2, 60.6, 63.1, 65.5,
    68.0, 70.5, 73.1, 75.7, 78.3
  ),
  q1 = c(
    .39278, .35913, .32882, .30128, .27608,
    .25289, .23143, .21151, .19292, .17553,
    .15920, .14383, .12912, .11524, .10259,
    .09037, .07862, .06734, .05656, .04628,
    .03604, .02678, .01838, .01118, .00581
  ),
  q2 = c(
    .47403, .43699, .40291, .37136, .34202,
    .31459, .28888, .26469, .24187, .22027,
    .19980, .18036, .16099, .14247, .12579,
    .10972, .09415, .07941, .06547, .05229,
    .03979, .02908, .01960, .01173, .00594
  ),
  q3 = c(
    .51004, .47149, .43575, .40242, .37123,
    .34193, .31433, .28825, .26354, .24010,
    .21780, .19654, .17511, .15453, .13611,
    .11843, .10138, .08521, .06989, .05537,
    .04178, .03033, .02030, .01204, .00611
  ),
  q5 = c(
    .55103, .51077, .47312, .43777, .40449,
    .37306, .34330, .31507, .28823, .26267,
    .23827, .21496, .19119, .16826, .14795,
    .12854, .11001, .09233, .07545, .05934,
    .04440, .03202, .02124, .01248, .00627
  ),
  q10 = c(
    .58262, .54211, .50388, .46771, .43339,
    .40077, .36969, .34005, .31171, .28459,
    .25861, .23368, .20814, .18342, .16141,
    .14033, .12015, .10083, .08236, .06469,
    .04830, .03475, .02297, .01342, .00669
  ),
  q15 = c(
    .60468, .56415, .52565, .48900, .45406,
    .42068, .38875, .35815, .32881, .30062,
    .27353, .24744, .22061, .19460, .17142,
    .14914, .12778, .10729, .08765, .06883,
    .05143, .03698, .02442, .01425, .00708
  ),
  q20 = c(
    .63218, .59182, .55318, .51613, .48057,
    .44640, .41354, .38189, .35139, .32198,
    .29358, .26614, .23796, .21061, .18593,
    .16214, .13924, .11721, .09604, .07571,
    .05675, .04093, .02711, .01587, .00793
  )
)

lifedeath_levels_north <- data.frame(
  level = 1:25,
  e0 = c(
    18.8, 21.2, 23.6, 26.1, 28.6,
    31.0, 33.5, 36.0, 38.4, 40.9,
    43.4, 45.9, 48.3, 50.8, 53.2,
    55.7, 58.1, 60.6, 63.1, 65.7,
    68.2, 70.7, 73.3, 75.9, 78.6
  ),
  q1 = c(
    .34608, .31620, .28938, .26509, .24292,
    .22256, .20377, .18635, .17012, .15496,
    .14076, .12744, .11503, .10386, .09300,
    .08249, .07235, .06262, .05331, .04445,
    .03602, .02806, .02094, .01417, .00922
  ),
  q2 = c(
    .42815, .39420, .36313, .33451, .30800,
    .28333, .26026, .23863, .21829, .19909,
    .18095, .16379, .14701, .13158, .11670,
    .10215, .08837, .07540, .06317, .05167,
    .04086, .03067, .02247, .01493, .00958
  ),
  q3 = c(
    .47836, .44191, .40825, .37698, .34782,
    .32050, .29482, .27062, .24775, .22609,
    .20553, .18602, .16657, .14861, .13139,
    .11473, .09899, .08414, .07015, .05696,
    .04451, .03271, .02369, .01557, .00963
  ),
  q5 = c(
    .54455, .50481, .46772, .43297, .40030,
    .36949, .34037, .31277, .28658, .26166,
    .23793, .21533, .19235, .17113, .15096,
    .13182, .11367, .09645, .08012, .06462,
    .04989, .03575, .02555, .01655, .01037
  ),
  q10 = c(
    .60282, .56241, .52408, .48765, .45296,
    .41987, .38826, .35802, .32904, .30127,
    .27459, .24902, .22271, .19815, .17474,
    .15248, .13131, .11120, .09208, .07388,
    .05653, .03985, .02820, .01804, .01115
  ),
  q15 = c(
    .62896, .58864, .55011, .51323, .47790,
    .44400, .41144, .38016, .35005, .32106,
    .29313, .26625, .23862, .21264, .18779,
    .16407, .14147, .11991, .09938, .07980,
    .06110, .04309, .03041, .01939, .01195
  ),
  q20 = c(
    .65402, .61407, .57561, .53856, .50284,
    .46841, .43517, .40308, .37207, .34209,
    .31311, .28513, .25654, .22941, .20332,
    .17831, .15436, .13145, .10953, .08858,
    .06851, .04912, .03476, .02242, .01399
  )
)

lifedeath_levels_south <- data.frame(
  level = 1:25,
  e0 = c(
    19.9, 22.4, 24.8, 27.2, 29.6,
    32.0, 34.4, 36.8, 39.2, 41.5,
    43.9, 46.3, 48.7, 51.0, 53.4,
    55.8, 58.1, 60.0, 63.1, 65.6,
    68.0, 70.5, 73.0, 75.5, 78.0
  ),
  q1 = c(
    .32162, .29818, .27701, .25772, .24004,
    .22373, .20862, .19455, .18141, .16929,
    .15862, .14821, .13806, .12819, .11861,
    .10934, .10040, .09170, .08178, .07215,
    .06256, .05304, .04363, .03442, .02658
  ),
  q2 = c(
    .45276, .42161, .39279, .36595, .34087,
    .31731, .29514, .27420, .25436, .23589,
    .21822, .21022, .18485, .16908, .15393,
    .13935, .12527, .11154, .09667, .08342,
    .07077, .05873, .04732, .03662, .02783
  ),
  q3 = c(
    .51319, .47850, .44614, .41583, .38733,
    .36045, .33501, .31090, .28798, .26658,
    .24569, .22564, .20640, .18793, .17020,
    .15318, .13678, .12082, .10382, .08898,
    .07492, .06168, .04929, .03781, .02852
  ),
  q5 = c(
    .56335, .52570, .49042, .45723, .42589,
    .39624, .36811, .34136, .31588, .29205,
    .26849, .24592, .22430, .20358, .18371,
    .16466, .14639, .12871, .11011, .09401,
    .07880, .06451, .05123, .03902, .02924
  ),
  q10 = c(
    .60047, .56231, .52613, .49175, .45901,
    .42777, .39790, .36932, .34191, .31575,
    .28999, .26525, .24149, .21868, .19677,
    .17570, .15547, .13589, .11578, .09836,
    .08201, .06676, .05270, .03990, .02974
  ),
  q15 = c(
    .61768, .57949, .54309, .50835, .47512,
    .44328, .41274, .38341, .35520, .32819,
    .30150, .27583, .25112, .22734, .20446,
    .18246, .16127, .14073, .11978, .10158,
    .08452, .06864, .05403, .04077, .03028
  ),
  q20 = c(
    .64117, .60305, .56648, .53133, .49752,
    .46495, .43356, .40327, .37402, .34606,
    .31809, .29110, .26506, .23993, .21571,
    .19235, .16983, .14795, .12559, .10626,
    .08815, .07135, .05593, .04200, .03105
  )
)

lifedeath_levels_east <- data.frame(
  level = 1:25,
  e0 = c(
    18.7, 21.2, 23.7, 26.2, 28.7,
    31.2, 33.6, 36.1, 38.6, 41.0,
    43.5, 45.9, 48.3, 50.7, 53.1,
    55.6, 58.0, 60.4, 62.8, 65.2,
    67.6, 70.0, 72.5, 75.0, 77.6
  ),
  q1 = c(
    .46740, .42976, .39546, .36396, .33488,
    .30791, .28277, .25925, .23719, .21643,
    .19701, .17907, .16179, .14515, .12917,
    .11384, .09915, .08513, .07176, .05902,
    .04631, .03488, .02460, .01576, .00886
  ),
  q2 = c(
    .53686, .49761, .46107, .42691, .39482,
    .36459, .33602, .30895, .28324, .25878,
    .23499, .21282, .19150, .17103, .15140,
    .13259, .11447, .09705, .08040, .06465,
    .05005, .03720, .02589, .01638, .00910
  ),
  q3 = c(
    .56608, .52614, .48867, .45337, .42002,
    .38842, .35841, .32984, .30260, .27658,
    .25096, .22700, .20399, .18191, .16075,
    .14047, .12095, .10221, .08428, .06726,
    .05186, .03836, .02656, .01671, .00923
  ),
  q5 = c(
    .59837, .55768, .51917, .48263, .44787,
    .41476, .38315, .35293, .32400, .29626,
    .26861, .24267, .21779, .19393, .17107,
    .14918, .12820, .10811, .08888, .07047,
    .05415, .03987, .02745, .01715, .00941
  ),
  q10 = c(
    .62691, .58639, .54768, .51065, .47515,
    .44110, .40839, .37692, .34663, .31744,
    .28806, .26042, .23384, .20827, .18371,
    .16013, .13748, .11573, .09487, .07485,
    .05741, .04210, .02883, .01792, .00976
  ),
  q15 = c(
    .64059, .60028, .56161, .52446, .48874,
    .45435, .42121, .38924, .35839, .32859,
    .29863, .27029, .24295, .21661, .19126,
    .16687, .14340, .12084, .09916, .07833,
    .06012, .04411, .03023, .01878, .01024
  ),
  q20 = c(
    .65908, .61918, .58067, .54349, .50757,
    .47282, .43921, .40667, .37514, .34460,
    .31400, .28478, .25651, .22918, .20279,
    .17733, .15277, .12911, .10632, .08437,
    .06490, .04780, .03291, .02056, .01129
  )
)

sbh_interpolation_row <- function(qxe_result, x, lifedeath_table, y){
  #define the qx (original x) and qy (y of interest to convert to) columns for filtering
  qx_col <- paste0('q', x)
  qy_col <- paste0('q', y)

  #pull relevant columns for interpolation from life death table
  qx_vals <- lifedeath_table[[qx_col]]
  qy_vals <- lifedeath_table[[qy_col]]

  #reverse orders as linear interpolation needs increasing x
  ord <- order(qx_vals)
  qx_vals <- qx_vals[ord]
  qy_vals <- qy_vals[ord]

  #perform linear interpolation to get the associated qy of interest
  qy_est <- stats::approx(
    x = qx_vals,
    y = qy_vals,
    xout = qxe_result,
    rule = 2
  )$y

  return(qy_est)
}


#' Conversion of q(x) estimates to a q(y) of interest
#'
#' This function takes q(x) estimates for a given x and uses Coale-Demeny life tables to convert it to a different q(y) of interest for that given reference year.
#' for more information on the process, refer to the demo vignette.
#'
#' @param sbh_result_table a table with previously estimated q(x)s and their reference years for all age groups/i for one specific region.
#' columns should be named region, q_i, i, x, and ref_year exactly. dataframe should be 7 rows long.
#'
#' @param family the family for the life death tables to be used (includes North, South, West, East - case sensitive)
#' @param y the under-y mortality rate of interest that all previous q(x)s will be mapped to.
#' this is NOT the same as i. for example, for u5mr you. input y = 5 and not y = 4. similarly, for under-10 you just input y = 10, not y = 5.
#'
#' @returns the sbh_result_table is returned with a new appended column representing the interpolated mortality rates.
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
#' sbh_interpolation_complete(
#'   sbh_result_table = sbh_results,
#'   family = "West",
#'   y = 5
#' )
#'
#' @export
sbh_interpolation_complete <- function(sbh_result_table, family, y) {
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
    lifedeath_table <- lifedeath_levels_west
  } else if (family == "North") {
    lifedeath_table <- lifedeath_levels_north
  } else if (family == "South") {
    lifedeath_table <- lifedeath_levels_south
  } else if (family == "East") {
    lifedeath_table <- lifedeath_levels_east
  }

  #for each q(x) in the summary table, compute the interpolated q(y)
  sbh_result_table |> dplyr::rowwise() |> dplyr::mutate(
    !!new_col := sbh_interpolation_row(
      qxe_result = q_i,
      x  = x,
      lifedeath_table = lifedeath_table,
      y = y
    )
  )
}
