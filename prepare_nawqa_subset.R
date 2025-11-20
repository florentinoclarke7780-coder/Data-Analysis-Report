## ------------------------------------------------------------
## NAWQA nutrients x land use: course-sized subset
## Purpose: Load original data, reproduce main joins, engineer
##          simple land-use summaries, aggregate to station level,
##          filter/sample a small dataset, and write CSV/RDS.
## ------------------------------------------------------------

## If needed, install packages (uncomment):
# install.packages(c("dplyr","readr","stringr","tidyr","here","ggplot2"))

## Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

## Optional helper for project-rooted paths
has_here <- requireNamespace("here", quietly = TRUE)
path_root <- function(...) {
  if (has_here) return(here::here(...))
  file.path(...)
}

message("[1/6] Loading raw data from CSV files in project root...")

## Reuse the authors' loading logic, adapting only paths.
## Original Rmd used base::read.csv with absolute paths. Here we use readr::read_csv
## on the repo-embedded CSVs with the same names.

## Nutrients (full NAWQA export)
NAWQA.NP <- readr::read_csv(
  file = path_root("NAWQA.NP.csv"),
  show_col_types = FALSE
)

## Station coordinates
nawqa.latlong <- readr::read_csv(
  file = path_root("nawqa.latlong.csv"),
  show_col_types = FALSE
)

## Watershed land use (StreamCat)
nawqa_landuse <- readr::read_csv(
  file = path_root("nawqa_landuse.csv"),
  show_col_types = FALSE
)

## After this section, we have three data frames in memory:
##   - NAWQA.NP       (nutrients)
##   - nawqa.latlong  (station coordinates)
##   - nawqa_landuse  (watershed land use metrics)

## ------------------------------------------------------------
## Reproduce the main join used by the authors
## (see nawqa_EAP_rmd_2020_02_26.Rmd)
## Steps mirrored:
## 1) Subset nutrient fields and rename to concise names (TN, TP, NO3, NH4, DIP)
## 2) Inner-join nutrients with coordinates by stationId/placeName/statePostalCode
## 3) Inner-join with land-use by latitude/longitude
## ------------------------------------------------------------

message("[2/6] Selecting and renaming key nutrient variables...")

## Handle possible differences in column names due to base::read.csv vs readr::read_csv
## In the CSV header, nitrate is "NO3+NO2_wf_00631"; in the Rmd examples it appears as
## "NO3.NO2_wf_00631". We coalesce both to a common NO3 variable.

pick <- function(df, name_candidates) {
  nm <- name_candidates[name_candidates %in% names(df)]
  if (length(nm) == 0) return(rep(NA_real_, nrow(df)))
  if (length(nm) == 1) return(df[[nm]])
  ## if multiple exist, take the first non-missing across them
  out <- df[[nm[1]]]
  if (length(nm) > 1) {
    for (k in nm[-1]) {
      out <- dplyr::coalesce(out, df[[k]])
    }
  }
  out
}

## Build a nutrient-only data frame, mirroring the Rmd approach
nutes_only <- tibble(
  stationId        = NAWQA.NP$stationId,
  statePostalCode  = NAWQA.NP$statePostalCode,
  placeName        = NAWQA.NP$placeName,
  resultDatetime   = NAWQA.NP$resultDatetime,
  TN               = NAWQA.NP$Totalnitrogen__62855,
  TN_code          = NAWQA.NP$Totalnitrogen__62855_RemarkCode,
  TP               = NAWQA.NP$Phosphorus_wu_00665,
  TP_code          = NAWQA.NP$Phosphorus_wu_00665_RemarkCode,
  NO3              = pick(NAWQA.NP, c("NO3.NO2_wf_00631","NO3+NO2_wf_00631")),
  NO3_code         = pick(NAWQA.NP, c("NO3.NO2_wf_00631_RemarkCode","NO3+NO2_wf_00631_RemarkCode")),
  NH4              = NAWQA.NP$Ammonia_wf_00608,
  NH4_code         = NAWQA.NP$Ammonia_wf_00608_RemarkCode,
  DIP              = NAWQA.NP$Orthophosphate__00671,
  DIP_code         = NAWQA.NP$Orthophosphate__00671_RemarkCode,
  siteVisit        = NAWQA.NP$siteVisitPurposeDescription
) %>%
  tidyr::drop_na(stationId, statePostalCode, placeName, TN, TP) # mirror na.omit on key fields

message("[3/6] Inner-joining nutrients with coordinates and land use...")

## Join nutrients + coordinates
nawqa_latlon_joined <- nutes_only %>%
  inner_join(
    nawqa.latlong,
    by = c("stationId","placeName","statePostalCode")
  )

## Join with land-use by latitude/longitude (as in the original Rmd)
nawqa_full <- nawqa_latlon_joined %>%
  inner_join(
    nawqa_landuse,
    by = c("latitude","longitude")
  )

## ------------------------------------------------------------
## Construct land-use summaries and a categorical land-use type
## ------------------------------------------------------------

message("[4/6] Engineering land-use summary variables and categories...")

## Use StreamCat watershed-percent variables (Ws). Fall back to Cat if Ws missing.
get_or_na <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NA_real_

nawqa_full <- nawqa_full %>%
  mutate(
    ## Summary percentages (use watershed-level variables available in StreamCat CSV)
    pAG  = PctCrop2006Ws + PctHay2006Ws,
    pURB = PctUrbLo2006Ws + PctUrbMd2006Ws + PctUrbHi2006Ws,
    pFOR = PctConif2006Ws + PctDecid2006Ws + PctMxFst2006Ws,

    ## Land use category heuristic (documented thresholds):
    ## - Agricultural: ag >= 50% and dominant
    ## - Urban:       urban >= 20% and dominant
    ## - Forest:      forest >= 50% and dominant
    ## - Otherwise:   Mixed
    LandUseType = dplyr::case_when(
      !is.na(pAG)  & pAG  >= 50 & pAG  >= pURB & pAG  >= pFOR ~ "Agricultural",
      !is.na(pURB) & pURB >= 20 & pURB >= pAG  & pURB >= pFOR ~ "Urban",
      !is.na(pFOR) & pFOR >= 50 & pFOR >= pAG  & pFOR >= pURB ~ "Forest",
      TRUE ~ "Mixed"
    )
  )

## ------------------------------------------------------------
## Reduce to one record per station (median of nutrients)
## ------------------------------------------------------------

message("[5/6] Aggregating to one record per station (medians)...")

## Helper to summarise if present, else NA
median_if_present <- function(df, var) {
  if (var %in% names(df)) median(.data[[var]], na.rm = TRUE) else NA_real_
}

station_level <- nawqa_full %>%
  group_by(
    stationId, statePostalCode, placeName,
    LandUseType, pAG, pURB, pFOR,
    latitude, longitude
  ) %>%
  summarise(
    TN  = median(TN,  na.rm = TRUE),
    TP  = median(TP,  na.rm = TRUE),
    NO3 = median(NO3, na.rm = TRUE),
    NH4 = median(NH4, na.rm = TRUE),
    .groups = "drop"
  )

## ------------------------------------------------------------
## Filter to a clean, small subset and sample up to ~120 rows
## ------------------------------------------------------------

message("[6/6] Filtering to clean subset and sampling up to 40 per land-use type...")

set.seed(202410)
subset_for_course <- station_level %>%
  filter(
    LandUseType %in% c("Agricultural", "Urban", "Forest"),
    !is.na(TN), TN > 0,
    !is.na(TP), TP > 0
  ) %>%
  group_by(LandUseType) %>%
  dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(40, nrow(.x)))) %>%
  ungroup()

## Final selection of columns for output
nawqa_course_data <- subset_for_course %>%
  select(
    stationId,
    statePostalCode,
    placeName,
    latitude,
    longitude,
    LandUseType,
    pAG,
    pURB,
    pFOR,
    TN,
    TP,
    NO3,
    NH4
  )

## Write outputs
out_dir <- path_root("data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

csv_path <- path_root("data","nawqa_landuse_nutrients_subset.csv")
rds_path <- path_root("data","nawqa_landuse_nutrients_subset.rds")

readr::write_csv(nawqa_course_data, csv_path)
try({ saveRDS(nawqa_course_data, rds_path) }, silent = TRUE)

## Summary table by land-use type
summary_by_landuse <- nawqa_course_data %>%
  group_by(LandUseType) %>%
  summarise(
    n = dplyr::n(),
    mean_TN = mean(TN, na.rm = TRUE),
    sd_TN   = sd(TN, na.rm = TRUE),
    mean_TP = mean(TP, na.rm = TRUE),
    sd_TP   = sd(TP, na.rm = TRUE),
    .groups = "drop"
  )

## Print a quick summary to console
message("\nSummary by LandUseType (n, mean/sd TN & TP):")
print(summary_by_landuse)

message(sprintf("\nWrote: %s", csv_path))
if (file.exists(rds_path)) message(sprintf("Also wrote: %s", rds_path))
