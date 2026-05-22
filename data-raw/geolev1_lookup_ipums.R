geolev1_lookup_ipums <- read.csv(
  "data-raw/geolev1_lookup_ipums.csv",
  stringsAsFactors = FALSE
)

#make the col numeric
geolev1_lookup_ipums$Geolev1 <- as.numeric(
  geolev1_lookup_ipums$Geolev1
)

usethis::use_data(
  geolev1_lookup_ipums,
  overwrite = TRUE
)
