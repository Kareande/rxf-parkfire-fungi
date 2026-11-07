# KAT ANDERSON ASV RICHNESS ANALYSIS WITH R
# Last updated 2026: 

#----------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------- PREPPING WORKSPACE -------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

## Install/load packages
#install.packages("tidyverse")
library(tidyverse)      # data manipulation + ggplot2

## Set working directory
# wd <- "C://Users//Kat//VS_Code//rxf-parkfire-fungi"
# wd <- "C://Users//mulle//VS_Code//rxf-parkfire-fungi"
setwd(wd)

## Load in data
frq_tbl <- read.delim("data/rarified_feat_freq_full.tsv", sep = "\t", check.names = FALSE, skip = 1) |>
  as.data.frame() |>
  rename_with(~ paste0("sample_", .x, recycle0 = TRUE)) |>        # rename columns to not begin w/ integer
  rename(asv = "sample_#OTU ID") |>                               # rename ASV column to ASV
  pivot_longer(-asv, values_to = "abundance") |>                  # make tidy with one row per ASV-Sample
  rename_all(tolower) |>                                          # make all columns lowercase
  mutate(site = substr(name, 8, nchar(name) - 2),                 # extract site from name (ex: 10 from sample_10B1)
         subsite = substr(name, 8, nchar(name) - 1)) |>           # extract subsite from name (ex: 10B sample_10B1)
  mutate(month = substr(name, nchar(name), nchar(name))) |>       # extract collection month
  mutate(month = str_replace_all(month, c("1" = "Jan", "7" = "July"))) |>                  # reformat months
  mutate(treatment = str_replace_all(site, c("10" = "NoRxF", "11" = "RxF", "12" = "NoRxF", # reformat treatments
                                             "1" = "RxF", "2" = "NoRxF", "3" = "RxF",
                                             "4" = "NoRxF", "5" = "NoRxF", "6" = "RxF",
                                             "7" = "NoRxF", "8" = "RxF", "9" = "RxF"))) |>
  mutate(pair = str_replace_all(site, c("10" = "PairY", "11" = "PairZ", "12" = "PairZ",    # designate pairs
                                        "1" = "PairT", "2" = "PairT", "3" = "PairV",
                                        "4" = "PairV", "5" = "PairW", "6" = "PairW",
                                        "7" = "PairX", "8" = "PairX", "9" = "PairY"))) |>
  mutate(tmntpair = factor(site,                                                           # create factors
                           levels = c("2", "1", "4", "3", "5", "6", "7", "8", "10", "9", "12", "11"),
                           labels = c("NoRxF Pair T", "RxF Pair T", "NoRxF Pair V", "RxF Pair V",
                                      "NoRxF Pair W", "RxF Pair W", "NoRxF Pair X", "RxF Pair X",
                                      "NoRxF Pair Y", "RxF Pair Y", "NoRxF Pair Z", "RxF Pair Z"))) |>
  mutate(treatment = factor(treatment),
         pair = factor(pair),
         site = factor(site, levels = c("2", "1", "4", "3", "5", "6",                      # RxF before NoRxF site
                                        "7", "8", "10", "9", "12", "11"))) |>
  select(asv, site, subsite, treatment, month, pair, tmntpair, abundance) |>                  # re-order desired columns
  filter(frq_tbl, abundance > 0)
frq_tbl


#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------- EXPLORE ASV RICHNESS -------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

## Number of universal ASV
freq_tbl_mnth <- frq_tbl |>
  mutate(site = paste(site, treatment, pair, sep = "_")) |>
  summarise(.by = c(asv, site, month),  # group data
            abundance = sum(abundance)) |>
  pivot_wider(names_from = site, values_from = abundance, values_fill = 0)

sum(apply(filter(freq_tbl_mnth, month == "Jan") != 0, 1, all)) #40 ASV at all sites, Jan
sum(apply(filter(freq_tbl_mnth, month == "July") != 0, 1, all)) #22 ASV at all sites, July

isTRUE(
  nrow(filter(freq_tbl_mnth, month == "July"))
  ==
    length(unique(filter(freq_tbl_mnth, month == "July")$asv))
)

freq_tbl_sites <- frq_tbl |>
  mutate(site = paste(site, treatment, pair, sep = "_")) |>
  summarise(.by = c(asv, site),  # group data
            abundance = sum(abundance)) |>
  pivot_wider(names_from = site, values_from = abundance, values_fill = 0)

sum(apply(freq_tbl_sites != 0, 1, all)) #61 ASV at all sites

freq_tbl_rxf <- frq_tbl |>
  filter(treatment == "RxF") |>
  mutate(site = paste(site, treatment, pair, sep = "_")) |>
  summarise(.by = c(asv, site, treatment),  # group data
            abundance = sum(abundance)) |>
  pivot_wider(names_from = site, values_from = abundance, values_fill = 0)

sum(apply(filter(freq_tbl_rxf, treatment == "RxF") != 0, 1, all)) #96 ASV at all RxF sites

isTRUE(
  nrow(filter(freq_tbl_rxf, treatment == "RxF"))
  ==
    length(unique(filter(freq_tbl_rxf, treatment == "RxF")$asv))
)

freq_tbl_norxf <- frq_tbl |>
  filter(treatment == "NoRxF") |>
  mutate(site = paste(site, treatment, pair, sep = "_")) |>
  summarise(.by = c(asv, site, treatment),  # group data
            abundance = sum(abundance)) |>
  pivot_wider(names_from = site, values_from = abundance, values_fill = 0)

sum(apply(filter(freq_tbl_norxf, treatment == "NoRxF") != 0, 1, all)) #72 ASV at all NoRxF sites

isTRUE(
  nrow(filter(freq_tbl_norxf, treatment == "NoRxF"))
  ==
    length(unique(filter(freq_tbl_norxf, treatment == "NoRxF")$asv))
)

## Number of unique ASV
freq_tbl_tmnt <- frq_tbl |>
  summarise(.by = c(asv, treatment),  # group data
            abundance = sum(abundance)) |>
  pivot_wider(names_from = treatment, values_from = abundance, values_fill = 0)

sum(apply(freq_tbl_tmnt != 0, 1, all)) #852 ASV present in both treatments
nrow(filter(freq_tbl_tmnt, RxF != 0 & NoRxF == 0)) #1020 ASV only exist in RxF
nrow(filter(freq_tbl_tmnt, RxF == 0 & NoRxF != 0)) #1040 ASV only exist in NoRxF

## ASV across paired sites
for (pairs in unique(frq_tbl$pair)){

  filt_freq <- frq_tbl |>
    filter(pair == pairs) |>
    mutate(pair = paste(treatment, pair, sep = "_")) |>
    summarise(.by = c(asv, pair),
              abundance = sum(abundance)) |>
    pivot_wider(names_from = pair, values_from = abundance, values_fill = 0)

  tot <- length(unique(filt_freq$asv))
  p_name <- names(filt_freq)
  print(paste("--- Now comparing pairs", p_name[2], "and", p_name[3],
              "which have", tot, "unique ASV ---"))

  p_all <- sum(apply(filt_freq != 0, 1, all))
  print(paste(p_all, "ASV present in both paired sites - ", (p_all/tot)*100, "%"))

  p_p1 <- nrow(filter(filt_freq, filt_freq[2] != 0 & filt_freq[3] == 0))
  p_p2 <- nrow(filter(filt_freq, filt_freq[2] == 0 & filt_freq[3] != 0))
  print(paste(p_p1, "ASV present in only", p_name[2], (p_p1/tot)*100, "%"))
  print(paste(p_p2, "ASV present in only", p_name[3], (p_p2/tot)*100, "%"))

}
rm(filt_freq)
