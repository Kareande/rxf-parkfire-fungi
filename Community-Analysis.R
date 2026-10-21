# KAT ANDERSON AVITI SEQUENCE ANALYSIS WITH R
# Updated 2026: 

#-------------------------------------------------------------#
#-------------------- PREPPING WORKSPACE ---------------------#
#-------------------------------------------------------------#

## Install/load packages
#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("viridis")
library(tidyverse)   # data manipulation + ggplot2
library(vegan)       # PERMANOVA (adonis2)
library(viridis)     # colour palettes
library(broom)       # for piping into tidy
library(ggtext)      # formatting plot labels

## Set working directory
wd = "C://Users//Kat//VS_Code//Bioinformatics"
setwd(wd)

## Load metadata
meta <- read.delim("metadata_noctl.tsv", sep = "\t", check.names = FALSE) |>
  group_by(Site)
meta$Treatment <- factor(meta$Treatment, levels = c("RxF", "NoRxF"))
meta <- meta[-1,]
head(meta)


#------------------------------------#
#---- GENERATE COLORS FOR PLOTS -----#
#------------------------------------#

## Set seed for colors (repeatability)
seed <- 82526
set.seed(seed) 

## Sample n colors from ColorBrewer
n <- 13*4
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- sort(col_vector)
n_colors <- sample(col_vector, n)

## Generate pallet of n colors
par(mfrow=c(1,1),bg='white')
pie(rep(1,n), col=n_colors, main=paste0('Color Pallet ',seed,'
(',n,' Colors)'))

## Save generated pallet
dev.print(device = png, filename = paste0('Images/ColorPallet_',seed,'_',n,'Col.png'), width=672, height=672)
dev.off()


#-------------------------------------------------------------#
#------------------ COMMUNITY - mikropml --------------------#
#-------------------------------------------------------------#

## Read frequency table
count_tbl <- read.delim("rarified_feat_freq_full.tsv", sep = "\t", check.names = FALSE, skip = 1) |>
  as.data.frame() |>
  rename(ASV = "#OTU ID") |>                             # rename ASV column to ASV
  pivot_longer(-ASV,
               names_to = "subsite",
               values_to = "count")
glimpse(count_tbl)

## Read taxonomy/function table
taxon_tbl <- read.delim("export-R/traits/fungal_traits.tsv", sep = "\t", check.names = FALSE) |>
  as.data.frame()
tnames <- names(taxon_tbl)
taxon_tbl <- taxon_tbl[-c(1), ]
taxon_tbl <- taxon_tbl |>
  select(tnames[1], tnames[2], tnames[4], tnames[5], tnames[6],
         tnames[7], tnames[8], tnames[9], tnames[10]) |>
  rename(ASV = "feature-id") |>
  mutate(taxa = str_extract(taxon, "(?<=\\g__)\\w+")) |> 
  mutate(taxa = if_else(is.na(taxa), str_extract(taxon, "([^;]*$)"), taxa)) |>
  mutate(taxon = taxa) |>
  select(-taxa)
glimpse(taxon_tbl)

## Read metadata
metadata <- read.delim("metadata/metadata_noctl.tsv", sep = "\t", check.names = FALSE) |>
  rename_all(tolower) |>
  rename(subsite = "sampleid") |>
  select(subsite, collectmonth, treatment, tmntpair)
metadata <- metadata[-c(1), ]
glimpse(metadata)

## Join count and taxon tables
abnd_tbl <- inner_join(count_tbl, taxon_tbl, by = "ASV") |>
  summarize(.by = c(subsite, taxon), count = sum(count)) |>
  group_by(subsite) |>
  mutate(rel_abund = count / sum(count)) |>
  ungroup() |>
  select(-count) |>
  inner_join(taxon_tbl) |>
  inner_join(metadata)
glimpse(abnd_tbl)

## What genera are significantly different between treatments?
sig_genera <- abnd_tbl |>
  nest(data = -taxon) |>                                   # all except taxon in a "data" column
  mutate(test = map(.x = data,                             # iterates function over columns
                    ~ wilcox.test(rel_abund ~ treatment,   # run wilcox test on each taxon row
                                  data = .x) |>
                      tidy())) |>                          # make nested test column tidy
  unnest(test)|>                                           # un-nest test data
  mutate(p.adjusted = p.adjust(p.value, method = "BH")) |> #adjust for multiple comparison false-positives
  filter(p.adjusted < 0.05) |>
  select(taxon, p.adjusted)
glimpse(sig_genera) #136 significant genera

sig_gen_plot <- abnd_tbl |>
  inner_join(sig_genera, by = "taxon") |>
  mutate(rel_abund = 100 * (rel_abund + 1/20000),
         taxon = str_replace(taxon, "(.*)", "*\\1*"),
#         taxon = str_repalce(taxon, "\\*(.*)_unclassified\\", "Unclassified<br>*\\1*"),
         treatment = factor(treatment, levels = c("NoRxF", "RxF"))) |>
  ggplot(aes(x = rel_abund, y = taxon, color = treatment, fill = treatment)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.3),
              shape = 12) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
               geom = "pointrange", position = position_dodge(width = 0.8),
               color = "black", show.legend = FALSE) +
  scale_x_log10() +
  labs(x = "Relative Abundance (%)", y = NULL) +
  scale_color_manual(NULL,
                     breaks = c(RxF, NoRxF),
                     values = c("gray", "red"),
                     labels = c("NoRxF", "RxF")) +
  theme_classic() +
  theme(axis.text.y = element_markdown())

sig_gen_plot
ggsave("C://Users//Kat//VS_Code//rxf-parkfire-fungi//figures//sig_genera.tiff", width =6, height=4)

## NOTE:
# "([^;]*$)", $ matches from end, * matches zero or more so ; does not HAVE to exist,
# [^.] matches one character outside of this (.) character class (https://regex101.com/r/fA3pI7/1/substitution)
#"(?<=\\g__)\\w+", (?<=\\.) matches anything after (.), \\w+ matches all word chars after
#
# Benjaminy Hotchbuerg less conservative/strict than Bonferonni; preferred for microbiome data

#-------------------------------------------------------------#
#----------------- DIFFERENTIAL ABUNDANCE --------------------#
#-------------------------------------------------------------#
## Install ANCOM-BC2 https://bioc.r-universe.dev/ANCOMBC/doc/manual.html
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")
library(ANCOMBC)

## Set up ASV table
#row names must be taxa and column names must be samples
asvs <- read.delim("export-R\\asv-table\\feat_table_taxon_full.tsv", skip = 1, sep = "\t", check.names = FALSE) |>
  as.data.frame() |>                                            # make it a dataframe
  mutate(taxonomy = gsub(" ", "", taxonomy)) |>
  select(ASV, taxonomy, everything())
nrow(asvs) #3998
length(unique(asvs$taxonomy)) #1059
head(asvs)
#                               ASV                                                                                                                   taxonomy        10B1        10C1        10D1        10E1        10G1        10H1
#1 be940f3d627f5e7aee26e5f93995658c                         k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Rhynchogastremaceae;g__Papiliotrema         272           0         171          90         287          15
#2 ea878a2bb3050c29de504f4c0ad4ed24    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;g__Cryptococcus;s__Cryptococcus_watticus        3634          32          17         755         616         279
#3 1fa535421f0269e28fe60ea966b9413c    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Filobasidiaceae;g__Naganishia;s__Naganishia_randhawae        1607          60          12        2225         357         481
#4 d49fcfae6b96a7c21466bb803c8a4565                                                                                                                   k__Fungi           0           0         244           0           0           2
#5 97079bf11a49c19a13468b7bee42a8b3 k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Piskurozymaceae;g__Solicoccozyma;s__Solicoccozyma_terrea        1399        7346          29         150          90          14
#6 cc8bf5cb112fac12a039a3561bd92bcc                                                                  k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales           0           0           0           7           0           0
asvs_ancom <- pivot_longer(asvs, c(-taxonomy,-ASV)) |> #remove the ASV column to only have taxon and sample freqs
  select(-ASV) |>                                       #remove ASV column
  rename(SampleID = "name") |>                          #rename column
  pivot_wider(names_from = SampleID, values_from = value, values_fn = sum) |>  #sum duplicates
  as.data.frame()
row.names(asvs_ancom) <- asvs_ancom$taxonomy
asvs_ancom <- as.data.frame.matrix(asvs_ancom[ , - 1])
nrow(asvs_ancom) #1059
length(names(asvs_ancom)) #150
head(asvs_ancom)[1:5,1:3]
#                                                                                                                   taxonomy        10B1        10C1
#                         k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Rhynchogastremaceae;g__Papiliotrema         272           0
#    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;g__Cryptococcus;s__Cryptococcus_watticus        4018          41
#    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Filobasidiaceae;g__Naganishia;s__Naganishia_randhawae        1614          60
#                                                                                                                   k__Fungi         764         529

## Set up metadata
#column names must be samples, must filter out samples not in feature table
metadata <- as.data.frame(read.delim("export-R\\asv-table\\metadata_noctl.tsv", sep = "\t", check.names = FALSE))
metadata <- metadata[-c(1), ]        #remove qiime data format row
row.names(metadata) <- metadata$SampleID
metadata <- metadata[,-1]
nrow(metadata) #216
head(metadata)[1:5,1:9]
#       Site CollectMonth Treatment TmntPair TmntPairMonth Subsample ...  Qubit Conc (ng/mL) Pooled Vol Sample_or_Control
#10A1   10          Jan     NoRxF    PairY NoRxFPairYJan         A ...              116.00        0.9      experimental
#10B1   10          Jan     NoRxF    PairY NoRxFPairYJan         B ...               21.70        0.6      experimental
#10C1   10          Jan     NoRxF    PairY NoRxFPairYJan         C ...              140.00        0.7      experimental
#10D1   10          Jan     NoRxF    PairY NoRxFPairYJan         D ...                9.28        0.5      experimental
#10E1   10          Jan     NoRxF    PairY NoRxFPairYJan         E ...              145.00        0.7      experimental

## Run differential abundance
#column names of metadata must match the sample names
set.seed(100126)
out = ancombc2(data = asvs_ancom, taxa_are_rows = TRUE, tax_level = "Genus",
                meta_data = metadata, fix_formula = "Treatment",
                rand_formula = NULL, #p_adj_method = "Site",
                pseudo_sens = TRUE,
                prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                group = "Site", struc_zero = FALSE, neg_lb = FALSE,
                alpha = 0.05, n_cl = 1, verbose = TRUE,
                global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
                em_control = list(tol = 1e-5, max_iter = 100, verbose = FALSE),
                #lme_control = lme4::lmerControl(),
                mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                #trend_control = list(contrast =
                #                        list(matrix(c(1, 0, -1, 1),
                #                                    nrow = 2,
                #                                    byrow = TRUE)),
                #                    node = list(2), B = 1)
                )


########
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
  "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))

res_prim <- out$res

m = out$res
datatable(m)
datatable(as.matrix(m))


df_rxf <- res_prim |>
    dplyr::select(taxon, ends_with("TreatmentRxF")) 
df_fig_TreatmentRxF <- df_rxf |> 
    dplyr::filter(diff_TreatmentRxF == "TRUE") |> 
    dplyr::arrange(desc(lfc_TreatmentRxF)) |>
    dplyr::mutate(direct = ifelse(lfc_TreatmentRxF > 0, "Positive LFC", "Negative LFC"),
                  color = ifelse(diff_robust_TreatmentRxF, "aquamarine3", "black")) |>
    mutate(taxon_s = str_extract(taxon, "(?<=\\g__)\\w+")) |>
    mutate(taxon_s = if_else(is.na(taxon_s), str_extract(taxon, "(?<=\\f__)\\w+"), taxon_s)) |>
    mutate(taxon_s = if_else(is.na(taxon_s), str_extract(taxon, "(?<=\\o__)\\w+"), taxon_s)) |>
    mutate(taxon_s = if_else(is.na(taxon_s), str_extract(taxon, "(?<=\\;c__)\\w+"), taxon_s))
df_fig_TreatmentRxF$taxon = factor(df_fig_TreatmentRxF$taxon, levels = df_fig_TreatmentRxF$taxon)
df_fig_TreatmentRxF$direct = factor(df_fig_TreatmentRxF$direct, 
                           levels = c("Positive LFC", "Negative LFC"))
df_fig_TreatmentRxF <- arrange(df_fig_TreatmentRxF, df_fig_TreatmentRxF$lfc_TreatmentRxF)

fig_TreatmentRxF = df_fig_TreatmentRxF |>
    ggplot(aes(y = lfc_TreatmentRxF, x = taxon_s, fill = direct)) + 
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lfc_TreatmentRxF - se_TreatmentRxF, 
                      ymax = lfc_TreatmentRxF + se_TreatmentRxF), 
                  width = 0.2,
                  position = position_dodge(0.05),
                  color = "black") + 
    labs(x = NULL,
         y = "Log fold change", 
         title = "Log fold changes as one unit increase of TreatmentRxF") + 
    scale_fill_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
           axis.text.x = element_text(angle = 60, hjust = 1))
fig_TreatmentRxF
ggsave(file = paste0("Images/ANCOM.jpeg"),
units = c("in"),
width = 8,
height = 5.5,
dpi = 300,
fig_TreatmentRxF)

df_bmi = res_prim %>%
    dplyr::select(taxon, contains("bmi")) 
df_fig_bmi1 = df_bmi %>%
    dplyr::filter(diff_bmilean == 1 | 
                    diff_bmioverweight == 1) %>%
    dplyr::mutate(lfc1 = ifelse(diff_bmioverweight == 1, 
                                round(lfc_bmioverweight, 2), 0),
                  lfc2 = ifelse(diff_bmilean == 1, 
                                round(lfc_bmilean, 2), 0)) %>%
    tidyr::pivot_longer(cols = lfc1:lfc2, 
                        names_to = "group", values_to = "value") %>%
    dplyr::arrange(taxon)

df_fig_bmi2 = df_bmi %>%
    dplyr::filter(diff_bmilean == 1 | 
                    diff_bmioverweight == 1) %>%
    dplyr::mutate(lfc1 = ifelse(diff_robust_bmioverweight, 
                                "aquamarine3", "black"),
                  lfc2 = ifelse(diff_robust_bmilean, 
                                "aquamarine3", "black")) %>%
    tidyr::pivot_longer(cols = lfc1:lfc2, 
                        names_to = "group", values_to = "color") %>%
    dplyr::arrange(taxon)

df_fig_bmi = df_fig_bmi1 %>%
    dplyr::left_join(df_fig_bmi2, by = c("taxon", "group"))

df_fig_bmi$group = recode(df_fig_bmi$group, 
                          `lfc1` = "Overweight - Obese",
                          `lfc2` = "Lean - Obese")
df_fig_bmi$group = factor(df_fig_bmi$group, 
                          levels = c("Overweight - Obese",
                                     "Lean - Obese"))
  
lo = floor(min(df_fig_bmi$value))
up = ceiling(max(df_fig_bmi$value))
mid = (lo + up)/2
fig_bmi = df_fig_bmi %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_bmi













asvs_6mo <- read.delim("export-R\\asv-table\\feat_table_taxon_6mo.tsv", skip = 1, sep = "\t", check.names = FALSE)
asvs_12mo <- read.delim("export-R\\asv-table\\feat_table_taxon_12mo.tsv", skip = 1, sep = "\t", check.names = FALSE)
nrow(asvs_6mo) #2181
nrow(asvs_12mo) #3134
nrow(asvs_6mo) + nrow(asvs_12mo) #5315


https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121
https://bioconductor.org/packages//release/bioc/html/biomformat.html
https://biom-format.org/documentation/biom_conversion.html

#where our results were driven by a pair of features that were only detectable with ASVs
https://forum.qiime2.org/t/differential-abundance-da-analysis-before-or-after-feature-table-collapse-to-a-specific-taxonomical-level/17682/3
https://journals.asm.org/doi/10.1128/msystems.00065-20



#-------------------------------------------------------------#
#-------------------- COMMUNITY TRAITS -----------------------#
#-------------------------------------------------------------#
## Load in ASV/Taxon frequency table
asvs <- read.delim("export-R\\asv-table\\feat_table_taxon_full.tsv", skip = 1, sep = "\t", check.names = FALSE) |>
  as.data.frame() |>                                            # make it a dataframe
  rename_with(~ paste0("sample_", .x, recycle0 = TRUE)) |>      # rename columns to not begin w/ integer
  rename(ASV = "sample_ASV") |>
  rename(taxonomy = "sample_taxonomy") |>
  mutate(taxonomy = gsub(" ", "", taxonomy)) |>
  select(ASV, taxonomy, everything())
nrow(asvs) #3998
head(asvs)
#                               ASV                                                                                                                   taxonomy  sample_10B1 sample_10C1 sample_10D1 sample_10E1 sample_10G1 sample_10H1
#1 be940f3d627f5e7aee26e5f93995658c                         k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Rhynchogastremaceae;g__Papiliotrema         272           0         171          90         287          15
#2 ea878a2bb3050c29de504f4c0ad4ed24    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;g__Cryptococcus;s__Cryptococcus_watticus        3634          32          17         755         616         279
#3 1fa535421f0269e28fe60ea966b9413c    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Filobasidiaceae;g__Naganishia;s__Naganishia_randhawae        1607          60          12        2225         357         481
#4 d49fcfae6b96a7c21466bb803c8a4565                                                                                                                   k__Fungi           0           0         244           0           0           2
#5 97079bf11a49c19a13468b7bee42a8b3 k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Piskurozymaceae;g__Solicoccozyma;s__Solicoccozyma_terrea        1399        7346          29         150          90          14
#6 cc8bf5cb112fac12a039a3561bd92bcc                                                                  k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales           0           0           0           7           0           0

## Load in Taxon/Trait table
taxtraits <- read.delim("export-R\\traits\\fungal_traits.tsv", sep = "\t", check.names = FALSE) |>
  as.data.frame() |> 
  rename(taxonomy = "taxon") |> 
  rename(ASV = "feature-id") |> 
  dplyr::mutate(across(everything(), ~na_if(trimws(.), "")))
taxtraits <- taxtraits[-c(1), ]
asvs$taxonomy <- stringr::str_trim(asvs$taxonomy)
nrow(taxtraits) #3998
head(taxtraits)
#                               ASV                                                                                                                   taxonomy        confidence fungal_traits_genus      primary_lifestyle
#2 be940f3d627f5e7aee26e5f93995658c                         k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Rhynchogastremaceae;g__Papiliotrema 0.999999999530847        Papiliotrema           mycoparasite
#3 ea878a2bb3050c29de504f4c0ad4ed24    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;g__Cryptococcus;s__Cryptococcus_watticus  0.99999902523924        Cryptococcus unspecified_saprotroph
#4 1fa535421f0269e28fe60ea966b9413c    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Filobasidiaceae;g__Naganishia;s__Naganishia_randhawae   0.9328056929248          Naganishia unspecified_saprotroph
#5 d49fcfae6b96a7c21466bb803c8a4565                                                                                                                   k__Fungi 0.999999999999995                <NA>                   <NA>
#6 97079bf11a49c19a13468b7bee42a8b3 k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Piskurozymaceae;g__Solicoccozyma;s__Solicoccozyma_terrea 0.999995085795457       Solicoccozyma        soil_saprotroph
#7 cc8bf5cb112fac12a039a3561bd92bcc                                                                  k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales 0.999999525894924                <NA>                   <NA>

## Join tables (frequency/taxon and trait)
traits <- left_join(asvs, taxtraits) |>
  select("ASV", "taxonomy", "confidence", everything()) 
names(traits)
head(traits)
nrow(traits)
write_delim(traits, "asv_traits.tsv", delim = "\t")
write_delim(traits, "asv_traits.csv", delim = ",")

## Subset for the Jan (6mo) samples
traits_6mo <- traits |>
    select(-ends_with("7")) |>                    #remove July (12mo) samples
    filter(if_any(everything(), ~ !is.na(.))) |>  #remove rows if everything is NA
    filter(if_any(ends_with("1"), ~ .x != 0))     #remove rows if the ASV frequency is 0 across all samples
names(traits_6mo)
nrow(traits_6mo) #2181
write_delim(traits_6mo, "asv_traits_6mo.tsv", delim = "\t")
write_delim(traits_6mo, "asv_traits_6mo.csv", delim = ",")

## Subset for the July (12mo) samples
traits_12mo <- traits |>
    select(-ends_with("1")) |>                    #remove Jan (6mo) samples
    filter(if_any(everything(), ~ !is.na(.))) |>  #remove rows if everything is NA
    filter(if_any(ends_with("7"), ~ .x != 0))     #remove rows if the ASV frequency is 0 across all samples
names(traits_12mo)
nrow(traits_12mo) #3134
write_delim(traits_12mo, "asv_traits_12mo.tsv", delim = "\t")
write_delim(traits_12mo, "asv_traits_12mo.csv", delim = ",")


https://cran.r-project.org/web/packages/stringr/vignettes/regular-expressions.html
https://www.rdocumentation.org/packages/DT/versions/0.34.0/topics/datatable
https://github.com/bioc/ANCOMBC
https://github.com/bioc/ANCOMBC/blob/devel/R/ancombc2.R
https://bioconductor.org/packages//release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
https://bioc.r-universe.dev/ANCOMBC/doc/manual.html
