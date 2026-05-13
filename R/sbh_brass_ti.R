#' Compute t(i) for calculated q(x)'s from the Brass (Trussell Variant) method using SBH information
#'
#' @param region_dataframe a dataframe for a single region/area with SBH information for 5.yr age bands & both sexes from 15-19 up to 45-49, and a corresponding index column (i) for the age groups... columns should be named: region, i, num_women, total_both, dead_both. all less region should be numeric. the dataframe should be 7 rows.
#' @param i_brass i corresponding to the q(x) of interest for the reference year. for example, i = 1 (ages 15-19) corresponds to x = 1, i = 2 (ages 15-19) corresponds to x = 2, etc. - this should be a single integer ranging from 1-7 (inclusive).
#' @param family string representing family of Coale-Demeny life tables to be used (options: East, West, North, South).
#' @param census_year single integer representing the year of the census being used
#'
#' @returns the reference year for a given i corresponding to the q(x) of interest
#' @export
sbh_brass_ti <- function(region_dataframe, i_brass, family, census_year) {
  stopifnot(
    is.data.frame(region_dataframe),
    all(c("region", "num_women", "total_both", "dead_both", "i") %in%
          names(region_dataframe)),
    nrow(region_dataframe) == 7,
    is.numeric(region_dataframe$i),
    all(region_dataframe$i == 1:7),
    is.numeric(region_dataframe$num_women),
    is.numeric(region_dataframe$total_both),
    is.numeric(region_dataframe$dead_both),
    length(unique(region_dataframe$region)) == 1,

    length(i_brass) == 1,
    is.numeric(i_brass),
    i_brass %% 1 == 0,
    i_brass %in% 1:7,

    length(family) == 1,
    family %in% c("West", "North", "East", "South"),

    length(census_year) == 1,
    is.numeric(census_year),
    census_year %% 1 == 0,
    census_year >= 1000,
    census_year <= 9999
  )

  a <- b <- c <- NULL
  if (family == "North") {
    if (i_brass == 1)      { a <- 1.0921; b <- 5.4732; c <- -1.9672}
    else if (i_brass == 2) { a <- 1.3207; b <- 5.3751; c <- 0.2133}
    else if (i_brass == 3) { a <- 1.5996; b <- 2.6268; c <- 4.3701}
    else if (i_brass == 4) { a <- 2.0779; b <- -1.7908; c <- 9.4126}
    else if (i_brass == 5) { a <- 2.7705; b <- -7.3403; c <- 14.9352}
    else if (i_brass == 6) { a <- 4.1520; b <- -12.2448; c <- 19.2349}
    else if (i_brass == 7) { a <- 6.9650; b <- -13.9160; c <- 19.9542}
  } else if (family == "West") {
    if (i_brass == 1)      { a <- 1.0970; b <- 5.5628; c <- -1.9956}
    else if (i_brass == 2) { a <- 1.3062; b <- 5.5677; c <- 0.2962}
    else if (i_brass == 3) { a <- 1.5305; b <- 2.5528; c <- 4.8962}
    else if (i_brass == 4) { a <- 1.9991; b <- -2.4261; c <- 10.4282}
    else if (i_brass == 5) { a <- 2.7632; b <- -8.4065; c <- 16.1787}
    else if (i_brass == 6) { a <- 4.3468; b <- -13.2436; c <- 20.1990}
    else if (i_brass == 7) { a <- 7.5242; b <- -14.2013; c <- 20.0162}

  } else if (family == "East") {
    if (i_brass == 1)      { a <- 1.0959; b <- 5.5864; c <- -1.9949}
    else if (i_brass == 2) { a <- 1.2921; b <- 5.5897; c <- 0.3631}
    else if (i_brass == 3) { a <- 1.5021; b <- 2.4692; c <- 5.0927}
    else if (i_brass == 4) { a <- 1.9347; b <- -2.6419; c <- 10.8533}
    else if (i_brass == 5) { a <- 2.6197; b <- -8.9693; c <- 17.0981}
    else if (i_brass == 6) { a <- 4.1317; b <- -14.3550; c <- 21.8247}
    else if (i_brass == 7) { a <- 7.3657; b <- -15.8083; c <- 22.3005}

  } else if (family == "South") {
    if (i_brass == 1)      { a <- 1.0900; b <- 5.4443; c <- -1.9721}
    else if (i_brass == 2) { a <- 1.3079; b <- 5.5568; c <- 0.2021}
    else if (i_brass == 3) { a <- 1.5173; b <- 2.6755; c <- 4.7471}
    else if (i_brass == 4) { a <- 1.9399; b <- -2.2739; c <- 10.3876}
    else if (i_brass == 5) { a <- 2.6157; b <- -8.4819; c <- 16.5153}
    else if (i_brass == 6) { a <- 4.0794; b <- -13.8308; c <- 21.1866}
    else if (i_brass == 7) { a <- 7.1796; b <- -15.3880; c <- 21.7892}
  }

  if (is.null(a) || is.null(b) || is.null(c)) {
    stop("Missing coefficients: check family or i_brass input.")
  }

  #calculate parity for i = 1,2,3
  p1 <- region_dataframe %>% filter(i == 1) %>% mutate(p1 = total_both / num_women) %>% pull(p1)
  p2 <- region_dataframe %>% filter(i == 2) %>% mutate(p2 = total_both / num_women) %>% pull(p2)
  p3 <- region_dataframe %>% filter(i == 3) %>% mutate(p3 = total_both / num_women) %>% pull(p3)

  #get components sorted
  pi1_2 <- p1/p2
  pi2_3 <- p2/p3

  ti <- (a + b*pi1_2 + c*pi2_3)
  ref_year <- census_year - ti
  return(ref_year)
}
