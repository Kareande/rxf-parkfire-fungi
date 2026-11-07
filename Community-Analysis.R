# KAT ANDERSON AVITI SEQUENCE ANALYSIS WITH R
# Updated 2026: 

#----------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------- PREPPING WORKSPACE -------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

## Load packages
library(tidyverse)   # data manipulation + ggplot2
library(vegan)       # PERMANOVA (adonis2)
library(viridis)     # colour palettes
library(broom)       # for piping into tidy
library(ggtext)      # formatting plot labels
library(rstatix)     # stats package, anova_test
library(nortest)     # anderson-darling test
library(car)         # leverne test

## Set working directory
# wd <- "C://Users//Kat//VS_Code//rxf-parkfire-fungi"
# wd <- "C://Users//mulle//VS_Code//rxf-parkfire-fungi"
# setwd(wd)


#----------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------- COMPOSITE DATASET --------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

## Load in metadata
metadata <- read.delim("data/metadata_noctl.tsv", sep = "\t", check.names = FALSE) |>
  rename_all(tolower) |>                              # make all columns lowercase
  rename(subsite = "sampleid", pair = "tmntpair") |>  # rename sampleid column
  rename(month = "collectmonth") |>                   # rename sampleid column
  select(site, subsite, month, treatment, pair)             # select columns
metadata <- metadata[-c(1), ]                         # remove Qiime2 formatting column
glimpse(metadata) #subsite, metadata... (one row per subsite: 216)

## Create ASV abundance (count/observations) table from frequency table
abund_tbl <- read.delim("data/rarified_feat_freq_full.tsv", sep = "\t", check.names = FALSE, skip = 1) |>
  as.data.frame() |>
  rename(asv = "#OTU ID") |>                             # rename ASV column to asv
  pivot_longer(-asv,                                     # give subsites/samples their own columns
               names_to = "subsite",
               values_to = "abundance")
glimpse(abund_tbl) #asv, subsite, abundance (one row per asv*subsite: 2912*115=334,880)

## Create table of taxa/function assignments for ASVs
taxon_code <- "g"                                                       # taxon of interest (k, p, c, o, f, g, or s)
taxa_tbl <- read.delim("data/fungal_traits.tsv", sep = "\t", check.names = FALSE) |>
  as.data.frame() |>
  select("feature-id":decay_substrate_template, -confidence) |>         # select columns
  rename(asv = "feature-id", taxonomy = "taxon") |>                     # rename columns
  mutate(taxa = str_extract(taxonomy,
                            paste0("(?<=", taxon_code, "__)\\w+"))) |>  # extract taxon from taxonomy column
  mutate(taxa = if_else(is.na(taxa),
                        str_extract(taxonomy, "([^;]*$)"),
                        taxa)) |>                                       # if no taxon, extract next taxon
  select(asv, taxa, everything(), -taxonomy, -fungal_traits_genus)     # remove taxonomy column
taxa_tbl <- taxa_tbl[-c(1), ]                                           # remove Qiime2 formatting column
glimpse(taxa_tbl) #asv, taxa, functions... (one row per unique ASV: 3,998)
#unique(taxa_tbl$taxa)

composite <- abund_tbl |>
  left_join(metadata, by = "subsite") |>        # add metadata to rarefied asv, which has 115/216 subsites (334,880 × 7)
  left_join(taxa_tbl, by = "asv") |>            # add taxa/functional data to dataset (334,880 × 14)
  summarize(.by = c(-asv, -abundance),
            abundance = sum(abundance)) |>      # calculate abundance of taxa per subsite (68,195 × 13)
  mutate(subsite = str_extract(subsite, "\\d+\\w"),     # extract through all digits and one word char
         site = str_extract(subsite, "\\d+")) |> # extract through all digits
  select(taxa, site, subsite,
         treatment, pair, month,
         abundance, everything())            # re-order columns
glimpse(composite) #one row per unique taxa*subsite (68,195 × 13)

## Histogram
ggplot(composite, aes(x = abundance)) +
  geom_histogram()  # bad kertosis

## Determine taxon based on code above
if (taxon_code == "s") {
  taxon <- "Species"
} else if (taxon_code == "g") {
  taxon <- "Genus"
} else if (taxon_code == "f") {
  taxon <- "Family"
} else if (taxon_code == "o") {
  taxon <- "Order"
} else if (taxon_code == "c") {
  taxon <- "Class"
} else if (taxon_code == "p") {
  taxon <- "Phylum"
} else {
  taxon <- "Kingdom"
}

## Create function to determine plot parameters
# plt_params <- function(taxon_code, metric) {
#     if (taxon_code == "s") {
#     taxon <- "Species"
#   } else if (taxon_code == "g") {
#     taxon <- "Genus"
#   } else if (taxon_code == "f") {
#     taxon <- "Family"
#   } else if (taxon_code == "o") {
#     taxon <- "Order"
#   } else if (taxon_code == "c") {
#     taxon <- "Class"
#   } else if (taxon_code == "p") {
#     taxon <- "Phylum"
#   } else {
#     taxon <- "Kingdom"
#   }
#   plttit <- paste(taxon, metric, "Across Paired Treatments and Months \n")
#   xtit <- "\n Site Pair"
#   ytit <- paste(taxon, metric, "\n")
#   metric <- metric
# }

## Designate plot color schemes
ctheme <- c("#807d7d", "#ac2b2b")                                             # designate fill colors for treatments
ctheme2 <- c("#008080", "#CC5500")                                            # designate fill colors for months
ccap <- "#4b4a4a"                                                             # designate caption font color

#----------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------- TAXA RICHNESS ----------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

## Determine list of taxa per site
for (sites in unique(composite$site)){

  filt_freq <- composite |>
    filter(site == sites) |>
    summarise(.by = c(taxa, site, treatment, pair, month),  # group data
              abundance = sum(abundance)) |>
    pivot_wider(names_from = month, values_from = abundance)

  assign(paste0("s", sites), filt_freq)

  rm(filt_freq)
}

## Calculate/verify species richness (# of unique taxa per subsite)
rich_tbl <- composite |>
  summarise(.by = c(site, subsite, treatment, pair, month),  # group data
            richness = specnumber(abundance)) |>             # calculate taxa richness
  select(subsite, treatment, pair,                           # re-order columns
         month, richness, everything())
length(unique(filter(composite, subsite == "10B" & abundance > 0)$taxa)) #62
length(unique(filter(composite, subsite == "11F" & abundance > 0)$taxa)) #54

## Check for outliers (keep running until all identified/removed)
identify_outliers(rich_tbl, richness) # one outlier
#rich_tbl <- filter(rich_tbl, subsite != "3E")

## Test normality (distribution assumption)
rich <- rich_tbl$richness
shapiro.test(rich)            #p = 0.3058 (g)
ks.test(rich, "pnorm")        #p < 2.2e-16 (g)
ad.test(rich)                 #p = 0.3567 (g)
plot(density(rich))           #normalish (g)
ggplot(rich_tbl, aes(sample = richness)) + #almost linear (g)
  geom_qq() +
  geom_qq_line()

## Test homoscedasticity (variance assumption)
ggplot(rich_tbl, aes(site, richness, group = treatment, color = treatment)) +
  geom_jitter()
leveneTest(richness ~ site, data = rich_tbl) #p = 0.701 (g)
plot(rich)                      # variance fine (g)
boxplot(rich)                   # symmetrical (g)

## Histogram
metric <- "Richness"
plttit <- paste("Frequencies of", taxon, metric, "\n")   # designate plot title
xtit <- paste(metric)                             # designate x-title
ytit <- paste("Density \n")        # designate y-title
plt_type <- "hist"
breaks <- pretty(range(rich), n = nclass.FD(rich), min.n = 1) # Freedman-Diaconis algorithm
bwidth <- breaks[2] - breaks[1]                                   # calculate ideal binwidth
ggplot(rich_tbl) +
  geom_histogram(aes(x = richness, y = ..density..),              # plot histogram scaled by density
                 fill = "steelblue",
                 colour = "black",
                 binwidth = bwidth) +
  labs(title = plttit,
       x = xtit,
       y = ytit) +
  stat_function(fun = dnorm,                                      # graph normal distribution
                args = list(mean = mean(rich_tbl$richness),
                            sd = sd(rich_tbl$richness)))
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## Visualize richness per treatment pair over across timepoints, point range (alt to boxplot)
#inspired by https://www.youtube.com/watch?v=KoAJZYEOwSE&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=240
metric <- "Richness"
plttit <- paste(taxon, metric, "Across Paired Treatments and Months \n")      # designate plot title
xtit <- "\n Site Pair"                                                        # designate x-title
ytit <- paste(metric, "\n")
plt_type <- "pntrng"
rich_tbl |>
  mutate(month = factor(month,                                              # make month a factor
                        levels = c("Jan", "July"),                          # order the months
                        labels = c("Samples Collected 6mo Post-fire",       # renaming Jan for plots
                                   "Samples Collected 12mo Post-fire"))) |> # renaming July for plots
  mutate(pair = factor(pair,
                       levels = c("PairT", "PairV", "PairW",
                                  "PairX", "PairY", "PairZ"),
                       labels = c("Pair T", "Pair V", "Pair W",
                                  "Pair X", "Pair V", "Pair Z"))) |>
  ggplot(aes(x = pair, y = richness, color = treatment)) +
  facet_wrap(~ month, nrow = 1, scales = "free_x", strip.position = "bottom") +        # split plot by months
  stat_summary(fun.data = median_hilow,                                                # medians for points
               fun.args = list(conf.int = 0.5),                                        # confidence interval for range
               geom = "pointrange",                                                    # point-range plot
               position = position_dodge(width = 0.5)) +                               # plot treatments side-by-side
  labs(title = paste(plttit),                                       # plot title
       subtitle = NULL,                                              # plot subtitle
       x = NULL,                                                     # x-axis label
       y = paste0(ytit),                                             # y-axis label
       caption = NULL,
       tag = "Figure X \n",                                          # tag label
       color = NULL) +                                               # legend title
  scale_color_manual(name = NULL,                                    # legend title
                     breaks = c("NoRxF", "RxF"),                     # legend categories
                     labels = c("Wildfire Only", "RxF + Wildfire"),  # legend category titles
                     values = ctheme) +                              # line fill colors
  theme(text = element_text(size = 14),                              # Overall text size
        plot.title.position = "plot",                                # align title to plot
        plot.title = element_text(hjust = 0.5,                       # title justification
                                  size = rel(1.2),                   # title font size
                                  face = "bold"),                    # title font
        axis.title.x = element_text(face = "bold",                   # x-axis font face
                                    size = rel(1)),                  # x-axis font size
        axis.title.y = element_text(face = "bold",                   # y-axis font face
                                    size = rel(1)),                  # y-axis font size
        plot.caption.position = "plot",                              # align caption to plot
        plot.caption = element_markdown(hjust = 0,                   # caption justification (markdown for HTML)
                                        size = rel(0.6),             # caption font size
                                        color = ccap),               # caption font color
        legend.title = element_text(face = "bold",                   # legend face
                                    size = rel(1)),                  # legend font size
        strip.placement = "outside",                                 # place strip below axis labels
        strip.background = element_blank(),                          # remove facet strip fill
        strip.text = element_text(face = "bold",                     # strip font face
                                  colour = "black",                  # strip font color
                                  size = rel(1)))                    # strip font size
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## Designate plot labels
metric <- "Mean Richness"
plttit <- paste(metric, "of", taxon, "Across Paired Treatments and Months \n")      # designate plot title
xtit <- "\n Site Pair"                                                        # designate x-title
ytit <- paste(metric, "\n")

## Visualize mean richness per paired site across timepoints
plt_type <- "barplt"
rich_tbl |> summarise(.by = c(site, treatment, month, pair),       # group data
                      avg_rich = mean(richness),
                      se = sd(richness) / sqrt(length(richness))) |>
  mutate(month = factor(month,
                        levels = c("Jan", "July"),
                        labels = c("Samples Collected 6mo Post-fire",
                                   "Samples Collected 12mo Post-fire"))) |>
  mutate(pair = factor(pair,
                       levels = c("PairT", "PairV", "PairW",
                                  "PairX", "PairY", "PairZ"),
                       labels = c("T", "V", "W",
                                  "X", "Y", "Z"))) |>
  ggplot(aes(x = pair, y = avg_rich, fill = treatment)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(x = pair, ymin = avg_rich - se, ymax = avg_rich + se),
                colour = "black",
                alpha = 0.6,
                size = 0.6,
                position = position_dodge2(.9, padding = .6)) +
  facet_grid(. ~ month)  +
  scale_fill_manual(breaks = c("NoRxF", "RxF"),                      # color by month
                    values = ctheme,                                 # bar fill colors
                    labels = c("Wildfire Only", "RxF + Wildfire")) + # legend labels
  labs(title = paste0(plttit),                                       # plot title
       subtitle = NULL,                                              # plot subtitle
       x = paste0(xtit),                                             # x-axis label
       y = paste0(ytit),                                             # y-axis label
       caption = NULL,
       tag = "Figure X \n",                                          # tag label
       fill = NULL) +                                                # legend title
  theme(text = element_text(size = 17),                              # Overall text size
        plot.title.position = "plot",                                # align title to plot
        plot.title = element_text(hjust = 0.5,                       # title justification
                                  size = rel(1.2),                   # title font size
                                  face = "bold"),                    # title font
        axis.title.x = element_text(face = "bold",                   # x-axis font face
                                    size = rel(1)),                  # x-axis font size
        axis.title.y = element_text(face = "bold",                   # y-axis font face
                                    size = rel(1)),                  # y-axis font size
        plot.caption.position = "plot",                              # align caption to plot
        plot.caption = element_markdown(hjust = 0,                   # caption justification (markdown for HTML)
                                        size = rel(0.6),             # caption font size
                                        color = ccap),               # caption font color
        legend.title = element_text(face = "bold",                   # legend face
                                    size = rel(1)),                  # legend font size
        strip.background = element_rect(fill = ctheme[1],            # strip background color
                                        colour = "white"),           # strip border color
        strip.text = element_text(face = "bold",                     # strip font face
                                  size = rel(0.8),                   # strip font size
                                  colour = "white"))                 # strip font color
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## Visualize mean richness per treatment pair over across timepoints, dumbell plot
#and https://www.youtube.com/watch?v=qaksmQabMUI&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=236
plttit <- paste("Change in", metric, "of", taxon,
                "Over 6 Months Across Paired Treatments \n")      # designate plot title
plt_type <- "dmbll"
rich_tbl |> summarise(.by = c(treatment, month, pair),            # group data
                      avg.rich = mean(richness)) |>
  pivot_wider(names_from = c(month),                              # use months for column titles
              values_from = avg.rich) |>                          # widen for bump calcs (to offset label)
  rename(avg.rich_Jan = "Jan", avg.rich_July = "July") |>         # rename so pivot_longer works easier
  mutate(pair = factor(pair,
                       levels = c("PairT", "PairV", "PairW",
                                  "PairX", "PairY", "PairZ"),
                       labels = c("T", "V", "W",
                                  "X", "Y", "Z"))) |>             # rename so pivot_longer works easier
  mutate(bump_Jan = if_else(avg.rich_Jan > avg.rich_July,         # create bumps 2 to the left or right of points
                            avg.rich_Jan + 2,
                            avg.rich_Jan - 2),
         bump_July = if_else(avg.rich_Jan > avg.rich_July,        # create bumps to the right or left of points
                             avg.rich_July - 2,
                             avg.rich_July + 2)) |>
  pivot_longer(cols = c(avg.rich_Jan, avg.rich_July, bump_Jan, bump_July),
               names_to = c(".value", "month"),                     # collapse into abnd or bump columns
               names_sep = "_") |>
  ggplot(aes(x = avg.rich, y = pair, color = month)) +
  facet_wrap(~ treatment, nrow = 1, strip.position = "bottom") +    # split plots by treatment
  geom_line(color = "black", linewidth = 1.75, alpha = 0.35) +         # add semi-trans black line of specific size
  geom_point(size = 2.5) +                                          # add point of particular size
  geom_text(aes(label = round(avg.rich, digits = 0), x = bump)) +   # label point w/ value at bump position
  scale_color_manual(breaks = c("Jan", "July"),                     # color by month
                     values = ctheme2,                              # point fill colors
                     labels = c("January", "July")) +               # legend labels
  labs(title = paste0(plttit),                                      # plot title
       subtitle = NULL,                                             # plot subtitle
       x = paste0(ytit),                                            # x-axis label
       y = paste0(xtit),                                            # y-axis label
       caption = "<i>Each 'dumbell' represents one site, of which contained 9 samples.</i>",
       tag = "Figure X \n",                                         # tag label
       color = NULL) +                                              # legend title
  theme_classic() +                                                 # remove gray blocks from background
  theme(plot.title.position = "plot",                               # align title to plot
        plot.title = element_text(hjust = 0.5,                      # title justification
                                  face = "bold",                    # title font face
                                  family = NULL),                   # title font
        axis.title.x = element_text(face = "bold",                  # x-axis font face
                                    size = 14),                     # x-axis font size
        axis.title.y = element_text(face = "bold",                  # y-axis font face
                                    size = 14),                     # y-axis font size
        plot.caption.position = "plot",                             # align caption to plot
        plot.caption = element_markdown(hjust = 0,                  # caption justification (markdown for HTML)
                                        color = ccap),              # caption font color
        legend.title = element_text(face = "bold",                  # legend face
                                    size = 12),                     # legend font size
        strip.background = element_rect(fill = "white",             # strip background color
                                        colour = "black"),          # strip border color
        strip.text = element_text(face = "bold",                    # strip font face
                                  colour = "black",                 # strip font color
                                  size = rel(1.2)))                 # strip font size
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## PARAMETRIC Paired t-Test, NONPARAMETRIC Wilcoxon signed-rank test
#t.test: Is there is a difference between the mean differences in taxa richness for the treatments?
#wilcox: Is there a difference between the summed ranks of the differences in taxa richness for the treatments?
rich_ptable <- rich_tbl |>                                     # use richness table
  summarise(.by = c(treatment, pair),                          # group data by treatment and pair
            richness = sum(richness)) |>                       # count number of unique taxa for each treatment*pair
  pivot_wider(names_from = treatment, values_from = richness)  # seperate columns for treatment richness
rich_ptable
t.test(rich_ptable$NoRxF, rich_ptable$RxF, paired = TRUE)      #null accepted (g)
wilcox.test(rich_ptable$NoRxF, rich_ptable$RxF, paired = TRUE) #null accepted (g)

## PARAMETRIC Paired t-Test, NONPARAMETRIC Wilcoxon signed-rank test
#t.test: Is there is a difference between the mean differences in taxa richness for the months?
#wilcox: Is there a difference between the summed ranks of the differences in taxa richness for the months?
rich_ptable2 <- rich_tbl |>                               # use richness table
  summarise(.by = c(month, pair),                         # group data by month and pair
            richness = sum(richness)) |>                  # count number of unique taxa for each month*pair
  pivot_wider(names_from = month, values_from = richness) # seperate columns for month richness
rich_ptable2
t.test(rich_ptable2$Jan, rich_ptable2$July, paired = TRUE)      #null accepted (g)
wilcox.test(rich_ptable2$Jan, rich_ptable2$July, paired = TRUE) #null accepted (g)

## PARAMETRIC ANOVA repeated measures, NONPARAMETRIC Friedman Test
#no benefit in choosing t.test vs ANOVA
#anova: Is there a difference between the mean taxa richness for the treatments?
#friedman: Is there a difference between the treatment rank sums?
rich_gtable <- rich_tbl |>                          # use richness table
  summarise(.by = c(treatment, month, pair),        # group data by treatment and month
            richness = mean(richness))              # count number of unique taxa for each treatment*month
print(rich_gtable, n = 24)
aov_res <- anova_test(data = rich_gtable,           # fit a model, taxa richness in response to treatment
                      dv = richness,                # dependant variable is count_asv
                      wid = pair,                   # subject identified is site
                      within = c(month, treatment))
get_anova_table(aov_res)                            #null accepted (g)
rich_gtbl <- rich_gtable |>                         # convert to matrix for Friedman Test
  pivot_wider(names_from = c(treatment, month),     # seperate columns for treatment*month richness
              values_from = richness)
rich_gmtrx <- as.matrix(rich_gtbl[, 2:5])           # convert to matrix for Friedman Test
rich_gmtrx
friedman.test(rich_gmtrx)                           #null accepted (g)


#----------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------- SHANNON DIVERSITY --------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

## Calculate shannon diversity
shan_tbl <- composite |>
  filter(abundance > 0) |>
  summarise(.by = c(site, subsite, treatment, pair, month),         # group data
            shannon = diversity(abundance, index = "shannon")) |>   # calculate shannon diversity
  select(site, subsite, treatment, pair,                            # re-order columns
         month, shannon, everything())
shan_tbl #one row per subsite (115 × 6)

## Check for outliers (keep running until all identified/removed)
identify_outliers(shan_tbl, shannon) #none

## Test normality (distribution assumption)
shan <- shan_tbl$shannon
shapiro.test(shan)            #null rejected (g)
ks.test(shan, "pnorm")        #null rejected (g)
ad.test(shan)                 #null rejected (g)
plot(density(shan))           #non-normal (g)
ggplot(shan_tbl, aes(sample = shannon)) + #not linear(g)
  geom_qq() +
  geom_qq_line()

## Test homoscedasticity (variance assumption)
ggplot(shan_tbl, aes(site, shannon, group = treatment, color = treatment)) +
  geom_jitter()
leveneTest(shannon ~ site, data = shan_tbl) #null accepted (g)
plot(shan)                      # variance poor (g)
boxplot(shan)                   # non-symmetrical (g)

## Histogram
metric <- "Shannon Diversity"
plttit <- paste("Frequencies of", taxon, metric, "\n")   # designate plot title
xtit <- paste(metric)                             # designate x-title
ytit <- paste("Density \n")        # designate y-title                          # designate y-title
plt_type <- "hist"
breaks <- pretty(range(shan), n = nclass.FD(shan), min.n = 1)     # Freedman-Diaconis algorithm
bwidth <- breaks[2] - breaks[1]                                   # calculate ideal binwidth
ggplot(shan_tbl) +
  geom_histogram(aes(x = shannon, y = ..density..),               # plot histogram scaled by density
                 fill = "steelblue",
                 colour = "black",
                 binwidth = bwidth) +
  labs(title = plttit,
       x = xtit,
       y = ytit) +
  stat_function(fun = dnorm,                                      # graph normal distribution
                args = list(mean = mean(shan_tbl$shannon),
                            sd = sd(shan_tbl$shannon)))
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## Visualize shannon diversity per treatment pair over across timepoints, point range (alt to boxplot)
#inspired by https://www.youtube.com/watch?v=KoAJZYEOwSE&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=240
metric <- "Shannon Diversity"
plttit <- paste(metric, "of", taxon,
                "Across Paired Treatments and Months \n")      # designate plot title
xtit <- "\n Site Pair"                                        # designate x-title
ytit <- paste(metric, "\n")                            # designate y-title
plt_type <- "pntrng"
shan_tbl |>
  mutate(month = factor(month,                                              # make month a factor
                        levels = c("Jan", "July"),                          # order the months
                        labels = c("Samples Collected 6mo Post-fire",       # renaming Jan for plots
                                   "Samples Collected 12mo Post-fire"))) |> # renaming July for plots
  mutate(pair = factor(pair,
                       levels = c("PairT", "PairV", "PairW",
                                  "PairX", "PairY", "PairZ"),
                       labels = c("Pair T", "Pair V", "Pair W",
                                  "Pair X", "Pair V", "Pair Z"))) |> # renaming July for plots
  ggplot(aes(x = pair, y = shannon, color = treatment)) +
  facet_wrap(~ month, nrow = 1, scales = "free_x", strip.position = "bottom") +        # split plot by months
  stat_summary(fun.data = median_hilow,                                                # medians for points
               fun.args = list(conf.int = 0.5),                                        # confidence interval for range
               geom = "pointrange",                                                    # point-range plot
               position = position_dodge(width = 0.5)) +                               # plot treatments side-by-side
  labs(title = paste(plttit),                                       # plot title
       subtitle = NULL,                                              # plot subtitle
       x = NULL,                                                     # x-axis label
       y = paste0(ytit),                                             # y-axis label
       caption = NULL,
       tag = "Figure X \n",                                          # tag label
       color = NULL) +                                               # legend title
  scale_color_manual(name = NULL,                                    # legend title
                     breaks = c("NoRxF", "RxF"),                     # legend categories
                     labels = c("Wildfire Only", "RxF + Wildfire"),  # legend category titles
                     values = ctheme) +                              # line fill colors
  theme(text = element_text(size = 14),                              # Overall text size
        plot.title.position = "plot",                                # align title to plot
        plot.title = element_text(hjust = 0.5,                       # title justification
                                  size = rel(1.2),                   # title font size
                                  face = "bold"),                    # title font
        axis.title.x = element_text(face = "bold",                   # x-axis font face
                                    size = rel(1)),                  # x-axis font size
        axis.title.y = element_text(face = "bold",                   # y-axis font face
                                    size = rel(1)),                  # y-axis font size
        plot.caption.position = "plot",                              # align caption to plot
        plot.caption = element_markdown(hjust = 0,                   # caption justification (markdown for HTML)
                                        size = rel(0.6),             # caption font size
                                        color = ccap),               # caption font color
        legend.title = element_text(face = "bold",                   # legend face
                                    size = rel(1)),                  # legend font size
        strip.placement = "outside",                                 # place strip below axis labels
        strip.background = element_blank(),                          # remove facet strip fill
        strip.text = element_text(face = "bold",                     # strip font face
                                  colour = "black",                  # strip font color
                                  size = rel(1)))                    # strip font size
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## Designate plot labels
metric <- "Mean Shannon Diversity"
plttit <- paste(metric, "of", taxon, "Across Paired Treatments and Months \n")      # designate plot title
xtit <- "\n Site Pair"                                                        # designate x-title
ytit <- paste(metric, "\n")

## Visualize mean richness per paired site across timepoints
plt_type <- "barplt"
shan_tbl |> summarise(.by = c(site, treatment, month, pair),       # group data
                     avg_shan = mean(shannon),
                     se = sd(shannon) / sqrt(length(shannon))) |>
  mutate(month = factor(month,
                        levels = c("Jan", "July"),
                        labels = c("Samples Collected 6mo Post-fire",
                                   "Samples Collected 12mo Post-fire"))) |>
  mutate(pair = factor(pair,
                       levels = c("PairT", "PairV", "PairW",
                                  "PairX", "PairY", "PairZ"),
                       labels = c("T", "V", "W",
                                  "X", "Y", "Z"))) |>
  ggplot(aes(x = pair, y = avg_shan, fill = treatment)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(x = pair, ymin = avg_shan - se, ymax = avg_shan + se),
                colour = "black",
                alpha = 0.6,
                size = 0.6,
                position = position_dodge2(.9, padding = .6)) +
  facet_grid(. ~ month)  +
  scale_fill_manual(breaks = c("NoRxF", "RxF"),                      # color by month
                    values = ctheme,                                 # bar fill colors
                    labels = c("Wildfire Only", "RxF + Wildfire")) + # legend labels
  labs(title = paste0(plttit),                                       # plot title
       subtitle = NULL,                                              # plot subtitle
       x = paste0(xtit),                                             # x-axis label
       y = paste0(ytit),                                             # y-axis label
       caption = NULL,
       tag = "Figure X \n",                                          # tag label
       fill = NULL) +                                                # legend title
  theme(text = element_text(size = 17),                              # Overall text size
        plot.title.position = "plot",                                # align title to plot
        plot.title = element_text(hjust = 0.5,                       # title justification
                                  size = rel(1.2),                   # title font size
                                  face = "bold"),                    # title font
        axis.title.x = element_text(face = "bold",                   # x-axis font face
                                    size = rel(1)),                  # x-axis font size
        axis.title.y = element_text(face = "bold",                   # y-axis font face
                                    size = rel(1)),                  # y-axis font size
        plot.caption.position = "plot",                              # align caption to plot
        plot.caption = element_markdown(hjust = 0,                   # caption justification (markdown for HTML)
                                        size = rel(0.6),             # caption font size
                                        color = ccap),               # caption font color
        legend.title = element_text(face = "bold",                   # legend face
                                    size = rel(1)),                  # legend font size
        strip.background = element_rect(fill = ctheme[1],            # strip background color
                                        colour = "white"),           # strip border color
        strip.text = element_text(face = "bold",                     # strip font face
                                  size = rel(0.8),                   # strip font size
                                  colour = "white"))                 # strip font color
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## Visualize mean shannon diversity per treatment pair over across timepoints, dumbell plot
#and https://www.youtube.com/watch?v=qaksmQabMUI&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=236
plttit <- paste("Change in", metric, "of", taxon,
                "Over 6 Months Across Paired Treatments \n")      # designate plot title
plt_type <- "dmbll"
shan_tbl |> summarise(.by = c(treatment, month, pair),       # group data
                      avg.shan = mean(shannon)) |>
  pivot_wider(names_from = c(month),                                # use months for column titles
              values_from = avg.shan) |>                            # widen for bump calcs (to offset label)
  rename(avg.shan_Jan = "Jan", avg.shan_July = "July") |>            # rename so pivot_longer works easier
  mutate(pair = factor(pair,
                       levels = c("PairT", "PairV", "PairW",
                                  "PairX", "PairY", "PairZ"),
                       labels = c("T", "V", "W",
                                  "X", "Y", "Z"))) |>           # rename so pivot_longer works easier
  mutate(bump_Jan = if_else(avg.shan_Jan > avg.shan_July,  # create bumps 2 to the left or right of points
                            avg.shan_Jan + .05,
                            avg.shan_Jan - .05),
         bump_July = if_else(avg.shan_Jan > avg.shan_July, # create bumps 2 to the right or left of points
                             avg.shan_July - .05,
                             avg.shan_July + .05)) |>
  pivot_longer(cols = c(avg.shan_Jan, avg.shan_July, bump_Jan, bump_July),
               names_to = c(".value", "month"),                     # collapse into abnd or bump columns
               names_sep = "_") |>
  ggplot(aes(x = avg.shan, y = pair, color = month)) +
  facet_wrap(~ treatment, nrow = 1, strip.position = "bottom") +    # split plots by treatment
  geom_line(color = "black", linewidth = 1.75, alpha = 0.35) +         # add semi-trans black line of specific size
  geom_point(size = 2.5) +                                          # add point of particular size
  geom_text(aes(label = round(avg.shan, digits = 1), x = bump)) +   # label point w/ value at bump position
  scale_color_manual(breaks = c("Jan", "July"),                     # color by month
                     values = ctheme2,                              # point fill colors
                     labels = c("January", "July")) +               # legend labels
  labs(title = paste0(plttit),                                      # plot title
       subtitle = NULL,                                             # plot subtitle
       x = paste0(ytit),                                            # x-axis label
       y = paste0(xtit),                                            # y-axis label
       caption = "<i>Each 'dumbell' represents one site, of which contained 9 samples.</i>",
       tag = "Figure X \n",                                         # tag label
       color = NULL) +                                              # legend title
  theme_classic() +                                                 # remove gray blocks from background
  theme(plot.title.position = "plot",                               # align title to plot
        plot.title = element_text(hjust = 0.5,                      # title justification
                                  face = "bold",                    # title font face
                                  family = NULL),                   # title font
        axis.title.x = element_text(face = "bold",                  # x-axis font face
                                    size = 14),                     # x-axis font size
        axis.title.y = element_text(face = "bold",                  # y-axis font face
                                    size = 14),                     # y-axis font size
        plot.caption.position = "plot",                             # align caption to plot
        plot.caption = element_markdown(hjust = 0,                  # caption justification (markdown for HTML)
                                        color = ccap),              # caption font color
        legend.title = element_text(face = "bold",                  # legend face
                                    size = 12),                     # legend font size
        strip.background = element_rect(fill = "white",             # strip background color
                                        colour = "black"),          # strip border color
        strip.text = element_text(face = "bold",                    # strip font face
                                  colour = "black",                 # strip font color
                                  size = rel(1.2)))                 # strip font size
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## NONPARAMETRIC Wilcoxon signed-rank test
#t.test: Is there is a difference between the mean differences in taxa richness for the treatments?
#wilcox: Is there a difference between the summed ranks of the differences in taxa richness for the treatments?
shan_ptable <- shan_tbl |>                                    # use richness table
  summarise(.by = c(treatment, pair),                         # group data by treatment and pair
            shannon = mean(shannon)) |>                       # count number of unique taxa for each treatment*pair
  pivot_wider(names_from = treatment, values_from = shannon)  # seperate columns for treatment richness
shan_ptable
wilcox.test(shan_ptable$NoRxF, shan_ptable$RxF, paired = TRUE) #null accepted (g)

## NONPARAMETRIC Wilcoxon signed-rank test
#wilcox: Is there a difference between the summed ranks of the differences in taxa richness for the months?
shan_ptable2 <- shan_tbl |>                               # use richness table
  summarise(.by = c(month, pair),                         # group data by month and pair
            shannon = mean(shannon)) |>                   # count number of unique taxa for each month*pair
  pivot_wider(names_from = month, values_from = shannon)  # seperate columns for month richness
shan_ptable2
wilcox.test(shan_ptable2$Jan, shan_ptable2$July, paired = TRUE) #null accepted (g)

## NONPARAMETRIC Friedman Test
#friedman: Is there a difference between the treatment rank sums?
shan_gtable <- shan_tbl |>                        # use richness table
  summarise(.by = c(treatment, month, pair),      # group data by treatment and month
            shannon = mean(shannon))              # count number of unique taxa for each treatment*month
print(shan_gtable, n = 24)
shan_gtbl <- shan_gtable |>                       # convert to matrix for Friedman Test
  pivot_wider(names_from = c(treatment, month),   # seperate columns for treatment*month richness
              values_from = shannon)
shan_gmtrx <- as.matrix(shan_gtbl[, 2:5])         # convert to matrix for Friedman Test
shan_gmtrx
friedman.test(shan_gmtrx)                         #null accepted (g)


#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------- ABUNDANCE VISUALIZATIONS ---------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

## Calculate relative abundance of each taxon by subsite
rabnd_tbl <- composite |>
  filter(abundance > 0) |>
  summarise(.by = c(taxa, site, subsite, treatment, month, pair),       # group data
            abundance = sum(abundance)) |>                              # mean rel abundance (%)
  mutate(.by = subsite, rel.abund = ((abundance / sum(abundance)) * 100)) |>
  select(-abundance) # (6,623 × 7)
rabnd_tbl

mrabnd_tbl <- rabnd_tbl |>
  summarise(.by = c(site, subsite, treatment, month, pair),       # group data
            mean.rabnd = mean(rel.abund))

## Identify outliers
ext <- identify_outliers(mrabnd_tbl, mean.rabnd) # (9 × 7)
flt <- mrabnd_tbl |> left_join(ext, by = join_by(taxa, site, subsite, treatment,
                                                month, pair, rel.abund)) |> # (6,623 × 9)
  arrange(is.extreme)
#flt <- flt[c(1:229, 1056:6623), ] |> select(-is.outlier, -is.extreme)

## Test normality (distribution assumption)
mrabnd <- mrabnd_tbl$mean.rabnd
shapiro.test(mrabnd)            #null rejected (g)
ks.test(mrabnd, "pnorm")        #null rejected (g)
ad.test(mrabnd)                 #null rejected (g)
plot(density(mrabnd))           #non-normal (g)
ggplot(mrabnd_tbl, aes(sample = mean.rabnd)) + #not linear(g)
  geom_qq() +
  geom_qq_line()

## Test homoscedasticity (variance assumption)
ggplot(mrabnd_tbl, aes(site, mean.rabnd, group = treatment, color = treatment)) +
  geom_jitter()
leveneTest(mean.rabnd ~ site, data = mrabnd_tbl) #null accepted (g)
plot(mrabnd)                      # variance poor (g)
boxplot(mrabnd)                   # non-symmetrical (g)

## Plot histogram
metric <- "Mean Relative Abundance"
plttit <- paste("Frequencies of", taxon, "Percent", metric, "\n")   # designate plot title
xtit <- paste(metric, "(%)")                             # designate x-title
ytit <- paste("Density \n")        # designate y-title                          # designate y-title
plt_type <- "hist"
breaks <- pretty(range(mrabnd), n = nclass.FD(mrabnd), min.n = 1) # Freedman-Diaconis algorithm
bwidth <- breaks[2] - breaks[1]                                   # calculate ideal binwidth
ggplot(mrabnd_tbl) +
  geom_histogram(aes(x = mean.rabnd, y = ..density..),              # plot histogram scaled by density
                 fill = "steelblue",
                 colour = "black",
                 binwidth = bwidth) +
  labs(title = plttit,
       x = xtit,
       y = ytit) +
  stat_function(fun = dnorm,                                      # graph normal distribution
                args = list(mean = mean(mrabnd),
                            sd = sd(mrabnd)))
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## Provide plot parameters
metric <- "Mean Relative Abundance"
plttit <- paste("Percent", metric, "of", taxon,
                "Across Paired Treatments and Months \n")      # designate plot title
xtit <- "\n Site Pair"                                         # designate x-title
ytit <- paste(metric, "(%)\n")                    # designate y-title

## Visualize abundance per treatment pair over across timepoints, point range (alt to boxplot)
#inspired by https://www.youtube.com/watch?v=KoAJZYEOwSE&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=240
plt_type <- "pntrng"
mrabnd_tbl |>
  mutate(month = factor(month,                                              # make month a factor
                        levels = c("Jan", "July"),                          # order the months
                        labels = c("Samples Collected 6mo Post-fire",       # renaming Jan for plots
                                   "Samples Collected 12mo Post-fire"))) |> # renaming July for plots
  ggplot(aes(x = pair, y = mean.rabnd, color = treatment)) +
  facet_wrap(~ month, nrow = 1, scales = "free_x", strip.position = "bottom") +        # split plot by months
  stat_summary(fun.data = median_hilow,                                                # medians for points
               fun.args = list(conf.int = 0.5),                                        # confidence interval for range
               geom = "pointrange",                                                    # point-range plot
               position = position_dodge(width = 0.5)) +                               # plot treatments side-by-side
  labs(title = paste(plttit),                                        # plot title
       subtitle = NULL,                                              # plot subtitle
       x = NULL,                                                     # x-axis label
       y = paste0(ytit),                                             # y-axis label
       caption = NULL,
       tag = "Figure X \n",                                          # tag label
       color = NULL) +                                               # legend title
  scale_color_manual(name = NULL,                                    # legend title
                     breaks = c("NoRxF", "RxF"),                     # legend categories
                     labels = c("Wildfire Only", "RxF + Wildfire"),  # legend category titles
                     values = ctheme) +                              # line fill colors
  theme(text = element_text(size = 14),                              # Overall text size
        plot.title.position = "plot",                                # align title to plot
        plot.title = element_text(hjust = 0.5,                       # title justification
                                  size = rel(1.2),                   # title font size
                                  face = "bold"),                    # title font
        axis.title.x = element_text(face = "bold",                   # x-axis font face
                                    size = rel(1)),                  # x-axis font size
        axis.title.y = element_text(face = "bold",                   # y-axis font face
                                    size = rel(1)),                  # y-axis font size
        plot.caption.position = "plot",                              # align caption to plot
        plot.caption = element_markdown(hjust = 0,                   # caption justification (markdown for HTML)
                                        size = rel(0.6),             # caption font size
                                        color = ccap),               # caption font color
        legend.title = element_text(face = "bold",                   # legend face
                                    size = rel(1)),                  # legend font size
        strip.placement = "outside",                                 # place strip below axis labels
        strip.background = element_blank(),                          # remove facet strip fill
        strip.text = element_text(face = "bold",                     # strip font face
                                  colour = "black",                  # strip font color
                                  size = rel(1)))                    # strip font size
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## Visualize mean richness per paired site across timepoints
plt_type <- "barplt"
rabnd_tbl |> summarise(.by = c(site, treatment, month, pair),       # group data
                       avg_rabnd = mean(rel.abund),
                       se = sd(rel.abund) / sqrt(length(rel.abund))) |>
  mutate(month = factor(month,
                        levels = c("Jan", "July"),
                        labels = c("Samples Collected 6mo Post-fire",
                                   "Samples Collected 12mo Post-fire"))) |>
  mutate(pair = factor(pair,
                       levels = c("PairT", "PairV", "PairW",
                                  "PairX", "PairY", "PairZ"),
                       labels = c("T", "V", "W",
                                  "X", "Y", "Z"))) |>
  ggplot(aes(x = pair, y = avg_rabnd, fill = treatment)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(x = pair, ymin = avg_rabnd - se, ymax = avg_rabnd + se),
                colour = "black",
                alpha = 0.6,
                size = 0.6,
                position = position_dodge2(.9, padding = .6)) +
  facet_grid(. ~ month)  +
  scale_fill_manual(breaks = c("NoRxF", "RxF"),                      # color by month
                    values = ctheme,                                 # bar fill colors
                    labels = c("Wildfire Only", "RxF + Wildfire")) + # legend labels
  labs(title = paste0(plttit),                                       # plot title
       subtitle = NULL,                                              # plot subtitle
       x = paste0(xtit),                                             # x-axis label
       y = paste0(ytit),                                             # y-axis label
       caption = NULL,
       tag = "Figure X \n",                                          # tag label
       fill = NULL) +                                                # legend title
  theme(text = element_text(size = 17),                              # Overall text size
        plot.title.position = "plot",                                # align title to plot
        plot.title = element_text(hjust = 0.5,                       # title justification
                                  size = rel(1.2),                   # title font size
                                  face = "bold"),                    # title font
        axis.title.x = element_text(face = "bold",                   # x-axis font face
                                    size = rel(1)),                  # x-axis font size
        axis.title.y = element_text(face = "bold",                   # y-axis font face
                                    size = rel(1)),                  # y-axis font size
        plot.caption.position = "plot",                              # align caption to plot
        plot.caption = element_markdown(hjust = 0,                   # caption justification (markdown for HTML)
                                        size = rel(0.6),             # caption font size
                                        color = ccap),               # caption font color
        legend.title = element_text(face = "bold",                   # legend face
                                    size = rel(1)),                  # legend font size
        strip.background = element_rect(fill = ctheme[1],            # strip background color
                                        colour = "white"),           # strip border color
        strip.text = element_text(face = "bold",                     # strip font face
                                  size = rel(0.8),                   # strip font size
                                  colour = "white"))                 # strip font color
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## Visualize mean abundance per treatment pair over across timepoints, dumbell plot
#and https://www.youtube.com/watch?v=qaksmQabMUI&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=236
plttit <- paste("Change in Percent", metric, "of", taxon,
                "Over 6 Months Across Paired Treatments \n")      # designate plot title
plt_type <- "dmbll"
rabnd_tbl |>
  summarise(.by = c(treatment, month, pair),             # group data by treatment, pair, and month
            avg.rel.abund = mean(rel.abund)) |>                    # mean relative abundance (%)
  pivot_wider(names_from = c(month),                                # use months for column titles
              values_from = avg.rel.abund) |>                            # widen for bump calcs (to offset label)
  rename(avg.rel.abund_Jan = "Jan", avg.rel.abund_July = "July") |>            # rename so pivot_longer works easier
  mutate(bump_Jan = if_else(avg.rel.abund_Jan > avg.rel.abund_July,  # create bumps 2 to the left or right of points
                            avg.rel.abund_Jan + .2,
                            avg.rel.abund_Jan - .2),
         bump_July = if_else(avg.rel.abund_Jan > avg.rel.abund_July, # create bumps 2 to the right or left of points
                             avg.rel.abund_July - .2,
                             avg.rel.abund_July + .2)) |>
  pivot_longer(cols = c(avg.rel.abund_Jan, avg.rel.abund_July, bump_Jan, bump_July),
               names_to = c(".value", "month"),                     # collapse into abnd or bump columns
               names_sep = "_") |>
  ggplot(aes(x = avg.rel.abund, y = pair, color = month)) +
  facet_wrap(~ treatment, nrow = 1, strip.position = "bottom") +    # split plots by treatment
  geom_line(color = "black", linewidth = 1.75, alpha = 0.35) +         # add semi-trans black line of specific size
  geom_point(size = 2.5) +                                          # add point of particular size
  geom_text(aes(label = round(avg.rel.abund, digits = 1), x = bump)) +   # label point w/ value at bump position
  scale_color_manual(breaks = c("Jan", "July"),                     # color by month
                     values = ctheme2,                              # point fill colors
                     labels = c("January", "July")) +               # legend labels
  labs(title = paste0(plttit),                                      # plot title
       subtitle = NULL,                                             # plot subtitle
       x = paste0(ytit),                                            # x-axis label
       y = paste0(xtit),                                            # y-axis label
       caption = "<i>Each 'dumbell' represents one site, of which contained 9 samples.</i>",
       tag = "Figure X \n",                                         # tag label
       color = NULL) +                                              # legend title
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
  theme_classic() +                                                 # remove gray blocks from background
  theme(plot.title.position = "plot",                               # align title to plot
        plot.title = element_text(hjust = 0.5,                      # title justification
                                  face = "bold",                    # title font face
                                  family = NULL),                   # title font
        axis.title.x = element_text(face = "bold",                  # x-axis font face
                                    size = 14),                     # x-axis font size
        axis.title.y = element_text(face = "bold",                  # y-axis font face
                                    size = 14),                     # y-axis font size
        plot.caption.position = "plot",                             # align caption to plot
        plot.caption = element_markdown(hjust = 0,                  # caption justification (markdown for HTML)
                                        color = ccap),              # caption font color
        legend.title = element_text(face = "bold",                  # legend face
                                    size = 12),                     # legend font size
        strip.background = element_rect(fill = "white",             # strip background color
                                        colour = "black"),          # strip border color
        strip.text = element_text(face = "bold",                    # strip font face
                                  colour = "black",                 # strip font color
                                  size = rel(1.2)))                 # strip font size
ggsave(paste0("figures/", taxon, "_", str_replace_all(metric, " ", "_"), "_", plt_type, ".png"))

## PARAMETRIC Paired t-Test, NONPARAMETRIC Wilcoxon signed-rank test
#t.test: Is there is a difference between the mean differences in taxa richness for the treatments?
#wilcox: Is there a difference between the summed ranks of the differences in taxa richness for the treatments?
rich_ptable <- rich_tbl |>                                     # use richness table
  summarise(.by = c(treatment, pair),                          # group data by treatment and pair
            richness = sum(richness)) |>                       # count number of unique taxa for each treatment*pair
  pivot_wider(names_from = treatment, values_from = richness)  # seperate columns for treatment richness
rich_ptable
t.test(rich_ptable$NoRxF, rich_ptable$RxF, paired = TRUE)      #null accepted (g)
wilcox.test(rich_ptable$NoRxF, rich_ptable$RxF, paired = TRUE) #null accepted (g)

## PARAMETRIC Paired t-Test, NONPARAMETRIC Wilcoxon signed-rank test
#t.test: Is there is a difference between the mean differences in taxa richness for the months?
#wilcox: Is there a difference between the summed ranks of the differences in taxa richness for the months?
rich_ptable2 <- rich_tbl |>                               # use richness table
  summarise(.by = c(month, pair),                         # group data by month and pair
            richness = sum(richness)) |>                  # count number of unique taxa for each month*pair
  pivot_wider(names_from = month, values_from = richness) # seperate columns for month richness
rich_ptable2
t.test(rich_ptable2$Jan, rich_ptable2$July, paired = TRUE)      #null accepted (g)
wilcox.test(rich_ptable2$Jan, rich_ptable2$July, paired = TRUE) #null accepted (g)

## PARAMETRIC ANOVA repeated measures, NONPARAMETRIC Friedman Test
#no benefit in choosing t.test vs ANOVA
#anova: Is there a difference between the mean taxa richness for the treatments?
#friedman: Is there a difference between the treatment rank sums?
rich_gtable <- rich_tbl |>                          # use richness table
  summarise(.by = c(treatment, month, pair),        # group data by treatment and month
            richness = mean(richness))              # count number of unique taxa for each treatment*month
print(rich_gtable, n = 24)
aov_res <- anova_test(data = rich_gtable,           # fit a model, taxa richness in response to treatment
                      dv = richness,                # dependant variable is count_asv
                      wid = pair,                   # subject identified is site
                      within = c(month, treatment))
get_anova_table(aov_res)                            #null accepted (g)
rich_gtbl <- rich_gtable |>                         # convert to matrix for Friedman Test
  pivot_wider(names_from = c(treatment, month),     # seperate columns for treatment*month richness
              values_from = richness)
rich_gmtrx <- as.matrix(rich_gtbl[, 2:5])           # convert to matrix for Friedman Test
rich_gmtrx
friedman.test(rich_gmtrx)                           #null accepted (g)











## Visualize mean asv counts per paired site across timepoints
pasv_plot <- avg_tbl_site |>
  mutate(month = factor(month,
                               levels = c("Jan", "July"),
                               labels = c("Samples Collected 6mo Post-fire",
                                          "Samples Collected 12mo Post-fire"))) |>
  ggplot(aes(x = pair, y = Avg_asv, fill = Treatment)) +
  geom_col(position = "dodge", alpha = 0.75) +
  facet_grid(. ~ CollectMonth)  +
  scale_fill_manual(values = c(col_vector[16], col_vector[22])) +
  labs(#title = "\n Mean ASV Abundance per Paired Site \n",
       x = "\n Site Pair",
       y = "Mean ASV Abundance\n") +
  theme(#plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = col_vector[16], colour = "white", size = 1),
        strip.text = element_text(face = "bold", colour = "white", size=rel(1.2)))
pasv_plot
ggsave("figures/Avg_ASV_Pair_Mo.png")


#----------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------- RELATIVE ABUNDANCE OF TAXA ----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
# inspired by https://www.youtube.com/watch?v=IoeSrLXeSos&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=251

## Calculate relative abundance of each taxon
rabnd_tbl <- composite |>
  mutate(.by = subsite, rel.abund = abundance / sum(abundance))  # calculate relative abundance
glimpse(rabnd_tbl)

## What taxon are significantly different? (all samples)
sig_taxa <- rabnd_tbl |>
  nest(data = -taxa) |>                                    # give each taxa attached data
  mutate(test = map(.x = data,                             # iterates function over columns
                    ~ anova_test(dv = rel.abund,           # dependant var
                                 wid = pair,               # identifier is site
                                 within = c(month, treatment),
                                 data = .x) |>             # run ANOVA on each taxon
                      tidy())) |>                          # make nested test column tidy
  unnest(test) |>                                          # un-nest test data
  mutate(p.adjusted = p.adjust(p.value, method = "BH")) |> # adjust for multiple comparison false-positives
  select(taxa, p.adjusted) |>                              # select only the taxon and adjusted p
  filter(p.adjusted < 0.05)
glimpse(sig_taxa)


## Visualize significant genera
sig_gen_plot <- rabnd_tbl |>
  inner_join(sig_taxa, by = "taxa") |>
  mutate(rel_abund = 100 * (rel.abund + 1/20000),
         taxa = str_replace(taxa, "(.*)", "*\\1*"),
#         taxon = str_repalce(taxa, "\\*(.*)_unclassified\\", "Unclassified<br>*\\1*"),
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
  scale_fill_manual(NULL,
                     breaks = c(RxF, NoRxF),
                     values = c("gray", "red"),
                     labels = c("NoRxF", "RxF")) +
  theme_classic() +
  theme(axis.text.y = element_markdown())

glimpse(sig_gen_plot)


sig_gen_plot
ggsave("C://Users//Kat//VS_Code//rxf-parkfire-fungi//figures//sig_genera.tiff", width =6, height=4)

## NOTE:
# "([^;]*$)", $ matches from end, * matches zero or more so ; does not HAVE to exist,
# [^.] matches one character outside of this (.) character class (https://regex101.com/r/fA3pI7/1/substitution)
#"(?<=\\g__)\\w+", (?<=\\.) matches anything after (.), \\w+ matches all word chars after
#
# Benjaminy Hotchbuerg less conservative/strict than Bonferonni; preferred for microbiome data

#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------- DIFFERENTIAL ABUNDANCE -----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
## Install ANCOM-BC2 https://bioc.r-universe.dev/ANCOMBC/doc/manual.html
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")
library(ANCOMBC)

## Set up ASV table
#row names must be taxa and column names must be samples
asvs <- read.delim("data/feat_table_taxon_full.tsv", skip = 1, sep = "\t", check.names = FALSE) |>
  as.data.frame() |>                                            # make it a dataframe
  mutate(taxonomy = gsub(" ", "", taxonomy)) |>
  select(asv, taxonomy, everything())
nrow(asvs) #3998
length(unique(asvs$taxonomy)) #1059
head(asvs)
#                               asv                                                                                                                   taxonomy        10B1        10C1        10D1        10E1        10G1        10H1
#1 be940f3d627f5e7aee26e5f93995658c                         k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Rhynchogastremaceae;g__Papiliotrema         272           0         171          90         287          15
#2 ea878a2bb3050c29de504f4c0ad4ed24    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;g__Cryptococcus;s__Cryptococcus_watticus        3634          32          17         755         616         279
#3 1fa535421f0269e28fe60ea966b9413c    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Filobasidiaceae;g__Naganishia;s__Naganishia_randhawae        1607          60          12        2225         357         481
#4 d49fcfae6b96a7c21466bb803c8a4565                                                                                                                   k__Fungi           0           0         244           0           0           2
#5 97079bf11a49c19a13468b7bee42a8b3 k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Piskurozymaceae;g__Solicoccozyma;s__Solicoccozyma_terrea        1399        7346          29         150          90          14
#6 cc8bf5cb112fac12a039a3561bd92bcc                                                                  k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales           0           0           0           7           0           0
asvs_ancom <- pivot_longer(asvs, c(-taxonomy,-asv)) |> #remove the asv column to only have taxon and sample freqs
  select(-asv) |>                                       #remove asv column
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
metadata <- as.data.frame(read.delim("data/metadata_noctl.tsv", sep = "\t", check.names = FALSE))
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













asvs_6mo <- read.delim("data/feat_table_taxon_6mo.tsv", skip = 1, sep = "\t", check.names = FALSE)
asvs_12mo <- read.delim("data/feat_table_taxon_12mo.tsv", skip = 1, sep = "\t", check.names = FALSE)
nrow(asvs_6mo) #2181
nrow(asvs_12mo) #3134
nrow(asvs_6mo) + nrow(asvs_12mo) #5315


#https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121
#https://bioconductor.org/packages//release/bioc/html/biomformat.html
#https://biom-format.org/documentation/biom_conversion.html

#where our results were driven by a pair of features that were only detectable with ASVs
#https://forum.qiime2.org/t/differential-abundance-da-analysis-before-or-after-feature-table-collapse-to-a-specific-taxonomical-level/17682/3
#https://journals.asm.org/doi/10.1128/msystems.00065-20



#----------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------- COMMUNITY TRAITS --------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
## Load in ASV/Taxon frequency table
asvs <- read.delim("data/feat_table_taxon_full.tsv", skip = 1, sep = "\t", check.names = FALSE) |>
  as.data.frame() |>                                            # make it a dataframe
  rename_with(~ paste0("sample_", .x, recycle0 = TRUE)) |>      # rename columns to not begin w/ integer
  rename(asv = "sample_asv") |>
  rename(taxonomy = "sample_taxonomy") |>
  mutate(taxonomy = gsub(" ", "", taxonomy)) |>
  select(asv, taxonomy, everything())
nrow(asvs) #3998
head(asvs)
#                               asv                                                                                                                   taxonomy  sample_10B1 sample_10C1 sample_10D1 sample_10E1 sample_10G1 sample_10H1
#1 be940f3d627f5e7aee26e5f93995658c                         k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Rhynchogastremaceae;g__Papiliotrema         272           0         171          90         287          15
#2 ea878a2bb3050c29de504f4c0ad4ed24    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;g__Cryptococcus;s__Cryptococcus_watticus        3634          32          17         755         616         279
#3 1fa535421f0269e28fe60ea966b9413c    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Filobasidiaceae;g__Naganishia;s__Naganishia_randhawae        1607          60          12        2225         357         481
#4 d49fcfae6b96a7c21466bb803c8a4565                                                                                                                   k__Fungi           0           0         244           0           0           2
#5 97079bf11a49c19a13468b7bee42a8b3 k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Piskurozymaceae;g__Solicoccozyma;s__Solicoccozyma_terrea        1399        7346          29         150          90          14
#6 cc8bf5cb112fac12a039a3561bd92bcc                                                                  k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales           0           0           0           7           0           0

## Load in Taxon/Trait table
taxtraits <- read.delim("data/fungal_traits.tsv", sep = "\t", check.names = FALSE) |>
  as.data.frame() |> 
  rename(taxonomy = "taxon") |> 
  rename(asv = "feature-id") |> 
  dplyr::mutate(across(everything(), ~na_if(trimws(.), "")))
taxtraits <- taxtraits[-c(1), ]
asvs$taxonomy <- stringr::str_trim(asvs$taxonomy)
nrow(taxtraits) #3998
head(taxtraits)
#                               asv                                                                                                                   taxonomy        confidence fungal_traits_genus      primary_lifestyle
#2 be940f3d627f5e7aee26e5f93995658c                         k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Rhynchogastremaceae;g__Papiliotrema 0.999999999530847        Papiliotrema           mycoparasite
#3 ea878a2bb3050c29de504f4c0ad4ed24    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;g__Cryptococcus;s__Cryptococcus_watticus  0.99999902523924        Cryptococcus unspecified_saprotroph
#4 1fa535421f0269e28fe60ea966b9413c    k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Filobasidiaceae;g__Naganishia;s__Naganishia_randhawae   0.9328056929248          Naganishia unspecified_saprotroph
#5 d49fcfae6b96a7c21466bb803c8a4565                                                                                                                   k__Fungi 0.999999999999995                <NA>                   <NA>
#6 97079bf11a49c19a13468b7bee42a8b3 k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Filobasidiales;f__Piskurozymaceae;g__Solicoccozyma;s__Solicoccozyma_terrea 0.999995085795457       Solicoccozyma        soil_saprotroph
#7 cc8bf5cb112fac12a039a3561bd92bcc                                                                  k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales 0.999999525894924                <NA>                   <NA>

## Join tables (frequency/taxon and trait)
traits <- left_join(asvs, taxtraits) |>
  select("asv", "taxonomy", "confidence", everything()) 
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


#https://cran.r-project.org/web/packages/stringr/vignettes/regular-expressions.html
#https://www.rdocumentation.org/packages/DT/versions/0.34.0/topics/datatable
#https://github.com/bioc/ANCOMBC
#https://github.com/bioc/ANCOMBC/blob/devel/R/ancombc2.R
#https://bioconductor.org/packages//release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
#https://bioc.r-universe.dev/ANCOMBC/doc/manual.html
