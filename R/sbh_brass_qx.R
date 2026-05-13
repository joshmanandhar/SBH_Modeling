#' Compute q(x) via Brass (Trussell Variant) method using SBH information
#'
#' @param district_dataframe a dataframe for a single region/area with SBH information for 5.yr age bands & both sexes from 15-19 up to 45-49, and a corresponding index column (i) for the age groups... columns should be named: region, i, num_women, total_both, dead_both. all less region should be numeric. the dataframe should be 7 rows.
#'
#' @param i_brass i corresponding to the q(x) being calculated. for example, i = 1 (ages 15-19) corresponds to x = 1, i = 2 (ages 15-19) corresponds to x = 2, etc. - this should be a single integer ranging from 1-7 (inclusive).
#' @param family string representing family of Coale-Demeny life tables to be used (options: East, West, North, South).
#'
#' @returns the region name from the dataframe, the computed q(x), and the multiplier k(x) from the associated formula
#' @export
sbh_brass_qx <- function(region_dataframe, i_brass, family) {
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
    family %in% c("West", "North", "East", "South")
  )
  #set coefficeints to be used and then map them to the proper ones based on the family and the i/age group
  a <- b <- c <- NULL
  if (family == "North") {
    if (i_brass == 1)      { a <- 1.1119; b <- -2.9287; c <- 0.8507}
    else if (i_brass == 2) { a <- 1.2390; b <- -0.6865; c <- -0.2745}
    else if (i_brass == 3) { a <- 1.1884; b <- 0.0421; c <- -0.5156}
    else if (i_brass == 4) { a <- 1.2046; b <- 0.3037; c <- -0.5656}
    else if (i_brass == 5) { a <- 1.2586; b <- 0.4236; c <- -0.5898}
    else if (i_brass == 6) { a <- 1.2240; b <- 0.4222; c <- -0.5456}
    else if (i_brass == 7) { a <- 1.1772; b <- 0.3486; c <- -0.4624}
  } else if (family == "West") {
    if (i_brass == 1)      { a <- 1.1415; b <- -2.7070; c <- 0.7663}
    else if (i_brass == 2) { a <- 1.2563; b <- -0.5381; c <- -0.2637}
    else if (i_brass == 3) { a <- 1.1851; b <- 0.0633; c <- -0.4177}
    else if (i_brass == 4) { a <- 1.1720; b <- 0.2341; c <- -0.4272}
    else if (i_brass == 5) { a <- 1.1865; b <- 0.3080; c <- -0.4452}
    else if (i_brass == 6) { a <- 1.1746; b <- 0.3314; c <- -0.4537}
    else if (i_brass == 7) { a <- 1.1639; b <- 0.3190; c <- -0.4435}
  } else if (family == "East") {
    if (i_brass == 1)      { a <- 1.1461; b <- -2.2536; c <- 0.6259}
    else if (i_brass == 2) { a <- 1.2231; b <- -0.4301; c <- -0.2245}
    else if (i_brass == 3) { a <- 1.1593; b <- 0.0581; c <- -0.3479}
    else if (i_brass == 4) { a <- 1.1404; b <- 0.1991; c <- -0.3487}
    else if (i_brass == 5) { a <- 1.1540; b <- 0.2511; c <- -0.3506}
    else if (i_brass == 6) { a <- 1.1336; b <- 0.2556; c <- -0.3428}
    else if (i_brass == 7) { a <- 1.1201; b <- 0.2362; c <- -0.3268}
  } else if (family == "South") {
    if (i_brass == 1)      { a <- 1.0819; b <- -3.0005; c <- 0.8689}
    else if (i_brass == 2) { a <- 1.2846; b <- -0.6181; c <- -0.3024}
    else if (i_brass == 3) { a <- 1.2223; b <- 0.0851; c <- -0.4704}
    else if (i_brass == 4) { a <- 1.1905; b <- 0.2631; c <- -0.4487}
    else if (i_brass == 5) { a <- 1.1911; b <- 0.3152; c <- -0.4291}
    else if (i_brass == 6) { a <- 1.1564; b <- 0.3017; c <- -0.3958}
    else if (i_brass == 7) { a <- 1.1307; b <- 0.2596; c <- -0.3538}
  }

  if (is.null(a) || is.null(b) || is.null(c)) {
    stop("Missing coefficients: check family or i_brass input.")
  }

  #calculate parity for i = 1,2,3 as neceesary for the formula
  p1 <- region_dataframe %>% filter(i == 1) %>% mutate(p1 = total_both / num_women) %>% pull(p1)
  p2 <- region_dataframe %>% filter(i == 2) %>% mutate(p2 = total_both / num_women) %>% pull(p2)
  p3 <- region_dataframe %>% filter(i == 3) %>% mutate(p3 = total_both / num_women) %>% pull(p3)

  #calculate d(i) for specified i
  di <- region_dataframe %>% filter(i == i_brass) %>% mutate(d_i = dead_both / total_both) %>% pull(d_i)

  #get components sorted
  pi1_2 <- p1/p2
  pi2_3 <- p2/p3

  #final formuala
  ki <- (a + b*pi1_2 + c*pi2_3)
  qi <- ki * di
  return(list('region' = region_dataframe$region[[1]],
              'ki' = ki,
              'qi' = qi))
}

