# KAT ANDERSON AVITI SEQUENCE ANALYSIS WITH R
# Last updated 2026: 

#----------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------- PREPPING WORKSPACE -------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

## Install/load packages
#install.packages("tidyverse")
#install.packages("RColorBrewer")
#install.packages("nortest")
#install.packages("vegan")
library(tidyverse)      # data manipulation + ggplot2
library(RColorBrewer)   # colour palettes
library(nortest)        # anderson-darling test
library(car)            # leverne test
library(vegan)          # microbio diversity
library(rstatix)        # stats package, anova_test
library(ggpubr)         # idk what this was for tbh
library(ggtext)         # element markdown to cinlude HTML
library(showtext)       # fonts
library(ggrepel)        # label text overlap

## Set working directory
wd <- "C://Users//Kat//VS_Code//rxf-parkfire-fungi"
#wd <- "C://Users//mulle//VS_Code//rxf-parkfire-fungi"
setwd(wd)

## Git repo
#git clone https://github.com/Kareande/rxf-parkfire-fungi.git
#git status
#git checkout septimelata
#git add .
#git reset HEAD
#git commit -m "comment message"    #-am skips stage step
#git commit --amend
#git push
#git pull origin septimelata


#----------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------- TIDY DATASET ----------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
## NOTE:
# I filtered out all samples with less than 500 reads

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
  select(asv, site, subsite, treatment, month, pair, tmntpair, abundance)                  # re-order desired columns
frq_tbl

## Histogram
ggplot(frq_tbl, aes(x = abundance)) +
  geom_histogram(binwidth = 20)  # yikes

## Remove any ASVs with zero counts
frq_tbl <- filter(frq_tbl, abundance > 0)

## Histogram
ggplot(frq_tbl, aes(x = abundance)) +
  geom_histogram(binwidth = 20)  # somewhat better


#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------- ABUNDANCE VISUALIZATIONS----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#plttit <- "\n Average ASV Abundance Across Paired Treatments and Months \n"
plttit <- "Average ASV Abundance Across Paired Treatments and Months \n" # designate plot title
xtit <- "\n Site Pair"                                                # designate x-title
ytit <- "Average ASV Abundance \n"                                    # designate y-title
ctheme <- c("#807d7d", "#ac2b2b")                                 # designate fill colors for treatments
ctheme2 <- c("#008080", "#CC5500")                                # designate fill colors for months
ccap <- "#4b4a4a"                                                   # designate caption font color

font_add_google(family = "patua-one", "Patua One")                    # download font
font_add_google(family = "roboto", "Roboto")                          # download font
showtext_auto(enable = TRUE)                                          # render the text
ftheme <- "patua-one"                                                 # designate font

## Visualize abundance per treatment pair over across timepoints, point range (alt to boxplot)
#inspired by https://www.youtube.com/watch?v=KoAJZYEOwSE&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=240
frq_tbl |> mutate(month = factor(month,                                               # make month a factor
                                 levels = c("Jan", "July"),                           # order the months
                                 labels = c("Samples Collected 6mo Post-fire",        # renaming Jan for plots
                                            "Samples Collected 12mo Post-fire"))) |>  # renaming July for plots
  ggplot(aes(x = pair, y = abundance, color = treatment)) +
  facet_wrap(~ month, nrow = 1, scales = "free_x", strip.position = "bottom") +       # split plot by months
  stat_summary(fun.data = median_hilow,                                               # medians for points
               fun.args = list(conf.int = 0.5),                                       # confidence interval for range
               geom = "pointrange",                                                   # point-range plot
               position = position_dodge(width = 0.5)) +                              # plot treatments side-by-side
  labs(title = paste0(plttit),                                       # plot title
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
ggsave("figures/abnd_pntrngplt.png")

## Visualize mean abundance per treatment pair over across timepoints, dumbell plot
#and https://www.youtube.com/watch?v=qaksmQabMUI&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=236
frq_tbl |> summarise(.by = c(treatment, month, pair),               # group data by treatment, pair, and month
                     avg.abnd = mean(abundance)) |>                 # average abundance
  pivot_wider(names_from = c(month),                                # use months for column titles
              values_from = avg.abnd) |>                            # widen for bump calcs (to offset label)
  rename(abnd_Jan = "Jan", abnd_July = "July") |>                   # rename so pivot_longer works easier
  mutate(bump_Jan = if_else(abnd_Jan > abnd_July,                   # create bumps 2 to the left or right of points
                            abnd_Jan + 2,
                            abnd_Jan - 2),
         bump_July = if_else(abnd_Jan > abnd_July,                  # create bumps 2 to the right or left of points
                             abnd_July - 2,
                             abnd_July + 2)) |>
  pivot_longer(cols = c(abnd_Jan, abnd_July, bump_Jan, bump_July),
               names_to = c(".value", "month"),                     # collapse into abnd or bump columns
               names_sep = "_") |>
  ggplot(aes(x = abnd, y = pair, color = month)) +
  facet_wrap(~ treatment, nrow = 1, strip.position = "bottom") +    # split plots by treatment
  geom_line(color = "black", linewidth = 1.75, alpha = 0.35) +         # add semi-trans black line of specific size
  geom_point(size = 2.5) +                                          # add point of particular size
  geom_text(aes(label = round(abnd, digits = 0), x = bump)) +       # label point w/ abundance value at bump position
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
ggsave("figures/avg_abnd_dmbllplt.png")

## Visualize mean abundance per treatment pair over across timepoints, slope plot
#and https://www.youtube.com/watch?v=z7Ci9NN3LUc&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=226
frq_tbl |> summarise(.by = c(treatment, month, pair),                 # group data by treatment, pair, and month
                     avg.abnd = mean(abundance)) |>
  mutate(treatment = factor(treatment,
                            levels = c("RxF", "NoRxF"),
                            labels = c("Prescibed Fire + Wildfire",
                                      "Wildfire Only"))) |>
  ggplot(aes(x = month, y = avg.abnd, group = pair, color = pair)) +
  facet_wrap(~ treatment, nrow = 1, strip.position = "bottom") +      # split plots by treatment
  geom_line() +
  labs(title = paste0(plttit),
       subtitle = NULL,
       caption = "<i>Each line represents one site and its change in mean abundance over 6 months.</i>",
       tag = "Figure X \n",
       x = NULL,
       y = paste0(ytit),
       color = "Treatment Paired Site") +
  theme(plot.title.position = "plot",                                 # align title to plot
        plot.title = element_text(face = "bold",                      # make title bold
                                  #family = ftheme,                    # change title font
                                  margin = margin(b=10)),             # change margin space (instead of \n)
        plot.caption = element_markdown(hjust = 0, color = ccap),     # move caption to be left-justified
        plot.caption.position = "plot",
        strip.placement = "outside",
        strip.background = element_blank())                           # align caption to plot
ggsave("figures/avg_abnd_slpplt.png")

## Visualize mean ASV counts per paired site
frq_tbl |>                                                           # avg observations
  summarise(.by = c(pair, treatment),
            avg_abnd = mean(abundance),
            se = sd(abundance) / sqrt(length(abundance))) |>
  ggplot(aes(x = pair, y = avg_abnd, fill = treatment)) +
  geom_col(position = "dodge", alpha = 0.9) +
  geom_errorbar(aes(x = pair, ymin = avg_abnd - se, ymax = avg_abnd + se),
                colour = "black",
                alpha = 0.6,
                size = 0.6,
                position = position_dodge2(.9, padding = .6)) +
  scale_fill_manual(values = ctheme) +
  labs(title = paste0(plttit),                                       # plot title
       subtitle = NULL,                                              # plot subtitle
       x = paste0(xtit),                                             # x-axis label
       y = paste0(ytit),                                             # y-axis label
       caption = NULL,
       tag = "Figure X \n",                                          # tag label
       fill = NULL) +                                         # legend title
  theme(text = element_text(size = 18),                              # Overall text size
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
                                    size = rel(1)))                  # legend font size
ggsave("figures/avg_abnd_barplt_fullstudy.png")

## Visualize mean ASV counts per paired site across timepoints
frq_tbl |> summarise(.by = c(site, treatment, month, pair, tmntpair),       # group data
                     avg_abnd = mean(abundance),
                     se = sd(abundance) / sqrt(length(abundance))) |>
  mutate(month = factor(month,
                        levels = c("Jan", "July"),
                        labels = c("Samples Collected 6mo Post-fire",
                                   "Samples Collected 12mo Post-fire"))) |>
  ggplot(aes(x = pair, y = avg_abnd, fill = treatment)) +
  geom_col(position = "dodge", alpha = 0.9) +
  geom_errorbar(aes(x = pair, ymin = avg_abnd - se, ymax = avg_abnd + se),
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
  theme(text = element_text(size = 18),                              # Overall text size
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
ggsave("figures/avg_abnd_barplt.png")


#----------------------------------------------------------------------------------------------------------------------#
#----------------------- DOES ASV RICHNESS DIFFER B/T TREATMENTS, MONTHS, OR TREATMENT*MONTHS? ------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#H0: The AVS counts do not differ across treatments
#H1: The AVS counts are different for soil treated and not treated with RxF

## Count unique ASVs per subsite/sample (richness)
rich_tbl <- frq_tbl |>                           # use freq table
  summarise(.by = c(site, subsite, treatment,    # group data
                    pair, tmntpair, month),
            count_asv = length(unique(asv)))     # number of unique ASV per subsite
rich_tbl

## Check for outliers (keep running until all identified/removed)
identify_outliers(rich_tbl, count_asv)
rich_tbl <- filter(rich_tbl, subsite != "3E")  #had to remove 3E, 12R, 3D
rich_tbl <- filter(rich_tbl, subsite != "12R") #had to remove 3E, 12R, 3D
rich_tbl <- filter(rich_tbl, subsite != "3D")  #had to remove 3E, 12R, 3D

## Test normality (distribution assumption)
shapiro.test(rich_tbl$count_asv)            #p = 0.5808, null accepted
ks.test(rich_tbl$count_asv, "pnorm")        #p < 2.2e-16, null rejected
ad.test(rich_tbl$count_asv)                 #p = 0.4196, null accepted
plot(density(rich_tbl$count_asv))           # plot density, normalish
ggplot(rich_tbl, aes(sample = count_asv)) + # plot Q-Q plot, almost linear
  geom_qq() +
  geom_qq_line()

## Test homoscedasticity (variance assumption)
ggplot(rich_tbl, aes(site, count_asv, group = treatment, color = treatment)) +
  geom_jitter()
leveneTest(count_asv ~ site, data = rich_tbl) #p = 0.7731, null accepted
plot(rich_tbl$count_asv)                      # plot variance, fine
boxplot(rich_tbl$count_asv)                   # plot spread, almost symmetrical

## PARAMETRIC Paired t-Test, NONPARAMETRIC Wilcoxon signed-rank test
#t.test: Is there is a difference between the mean differences in ASV counts for the treatments?
#wilcox: Is there a difference between the summed ranks of the differences in ASV counts for the treatments?
count_ptable <- rich_tbl |>                                      # use count table
  summarise(.by = c(treatment, pair),                            # group data by treatment and pair
            count_asv = sum(count_asv)) |>                       # count number of unique ASV for each treatment*pair
  pivot_wider(names_from = treatment, values_from = count_asv)   # seperate columns for treatment counts
count_ptable
t.test(count_ptable$NoRxF, count_ptable$RxF, paired = TRUE)      #p = 0.8132, null accepted
wilcox.test(count_ptable$NoRxF, count_ptable$RxF, paired = TRUE) #p = 0.6875, null accepted

## PARAMETRIC Paired t-Test, NONPARAMETRIC Wilcoxon signed-rank test
#t.test: Is there is a difference between the mean differences in ASV counts for the months?
#wilcox: Is there a difference between the summed ranks of the differences in ASV counts for the months?
count_ptable2 <- rich_tbl |>                               # use count table
  summarise(.by = c(month, pair),                          # group data by month and pair
            count_asv = sum(count_asv)) |>                 # count number of unique ASV for each month*pair
  pivot_wider(names_from = month, values_from = count_asv) # seperate columns for month counts
count_ptable2
t.test(count_ptable2$Jan, count_ptable2$July, paired = TRUE)      #p = 0.5118, null accepted
wilcox.test(count_ptable2$Jan, count_ptable2$July, paired = TRUE) #p = 0.5625, null accepted

## PARAMETRIC ANOVA repeated measures, NONPARAMETRIC Friedman Test
#no benefit in choosing t.test vs ANOVA
#anova: Is there a difference between the mean ASV counts for the treatments?
#friedman: Is there a difference between the treatment rank sums?
count_gtable <- rich_tbl |>                         # use count table
  summarise(.by = c(treatment, month, pair),        # group data by treatment and month
            count_asv = mean(count_asv))            # count number of unique ASV for each treatment*month
print(count_gtable, n = 24)
aov_res <- anova_test(data = count_gtable,          # fit a model, ASV richness in response to treatment
                      dv = count_asv,               # dependant variable is count_asv
                      wid = pair,                   # subject identified is site
                      within = c(month, treatment))
get_anova_table(aov_res)                            #p > 0.05, null accepted
count_gtbl <- count_gtable |>                       # convert to matrix for Friedman Test
  pivot_wider(names_from = c(treatment, month),     # seperate columns for treatment*month counts
              values_from = count_asv)
count_gmtrx <- as.matrix(count_gtbl[, 2:5])         # convert to matrix for Friedman Test
count_gmtrx
friedman.test(count_gmtrx)                          #p = 0.7055, null accepted

## NOTE:
# Data is normal enough to justify using parametric tests. Did both, but will only use parametric results.
# There is no support for the hypothesis that ASV richness is different between treatments, timepoints, nor
# treatments across timepoints.

## Test timepoint difference Roualdes' way
## Does the mean differences in ASV counts significantly differ across treatments (all sites) or pairs?
diff_table <- count_gtbl |>                                      # use Friedman Test table
  mutate(NoRxF_diff = NoRxF_July - NoRxF_Jan, pair,              # subtract Jan counts from July
         RxF_diff = RxF_July - RxF_Jan) |>
  select(pair, RxF_diff, NoRxF_diff) |>
  pivot_longer(c(RxF_diff, NoRxF_diff),
               names_to = "treatment",
               values_to = "diff_asv") |>
  mutate(treatment = substr(treatment, 1, nchar(treatment) - 5)) # trim the site from subsite
diff_table

## PARAMETRIC ANOVA, NONPARAMETRIC Kruskal-Wallis
#anova test sig diff between the means by analyzing variance
fit <- lm(diff_asv ~ treatment, data = diff_table)    # fit model, difference in response to tretament
anova(fit)                                            #p = 0.03106, null rejected

## Post-hoc test for above



## NOTE:
# Timepoint analysis can't be done by subsite because not all subsites have their 6 and
# 12mo sample comparisons for two reasons (1) the mislabeling of sites led to some sites
# not having their A-I match and instead being J-R, (2) only 157 samples out of 216 made
# it through processing.
#
# There is no support for the hypothesis that the richness is different across timepoints for either
# treatment (Wildfire only, or RxF + wildfire).
#
# Professor Roualdes initially grouped by site (summing the counts for each subsite),
# but then stated because they are their own "locations" and not some "smallest unit
# of measurement", that they don't need to be collapsed. He then said it's okay
# to remove the zeroes, and the ecologist stats guide supported this by saying that
# making zero observation counts important (that is, keeping them for what they are)
# conflates their meaning. After removing zeroes, the data was normal enough to use
# perimetric ANOVA. Still, no significant difference between diversity.
# Roualdes also says that "this analysis would be enough for me", but said he is not my
# comittee nor a peer reviewer so I could still look at the other metrics (Shannon diversity
# and eveness) in the same fashion.
# He also did not seperate out the timepoints, which at first troubled me.
# After thinking about it, for the question we were answering it seemed appropriate. If
# sampling over time is a way to deal with pseudoreplication (per the paper he gave us),
# then leaving it as one dataset seems more statistically sound than seprating them out and
# then having two, more pseudreplicated datasets than one. I still looked anyway, and it
# didn't make a difference.


#----------------------------------------------------------------------------------------------------------------------#
#---------------------- DOES SHANNON/SIMPSON DIFFER B/T TREATMENTS, MONTHS, OR TREATMENT*MONTHS? ----------------------#
#----------------------------------------------------------------------------------------------------------------------#

## Calculate diversity indexes
shn_tbl <- frq_tbl |>                                            # use ASV table
  summarise(.by = c(site, pair, subsite, treatment, month),      # group data
            shannon = diversity(abundance, index = "shannon"))
shn_tbl
smp_tbl <- frq_tbl |>                                            # use ASV table
  summarise(.by = c(site, pair, subsite, treatment, month),      # group data
            simpson = diversity(abundance, index = "simpson"))
smp_tbl

## Check for outliers (keep running until all identified/removed)
identify_outliers(shn_tbl, shannon)
identify_outliers(smp_tbl, simpson)
smp_tbl <- filter(smp_tbl, subsite != "6R") #had to remove 6R, 5G
smp_tbl <- filter(smp_tbl, subsite != "5G") #had to remove 6R, 5G

## Test normality (distribution assumption)
shapiro.test(shn_tbl$shannon)            #p = 0.0002574, null rejected
shapiro.test(smp_tbl$simpson)            #p = 5.899e-11, null rejected
ks.test(shn_tbl$shannon, "pnorm")        #p < 2.2e-16, null rejected
ks.test(smp_tbl$simpson, "pnorm")        #p < 2.2e-16, null rejected
ad.test(shn_tbl$shannon)                 #p = 4.176e-05, null rejected
ad.test(smp_tbl$simpson)                 #p < 2.2e-16, null rejected
plot(density(shn_tbl$shannon))           # plot density, left kertosis
plot(density(smp_tbl$simpson))           # plot density, left kertosis
ggplot(shn_tbl, aes(sample = shannon)) + # plot Q-Q plot, bad
  geom_qq() +
  geom_qq_line()
ggplot(smp_tbl, aes(sample = simpson)) + # plot Q-Q plot, even worse
  geom_qq() +
  geom_qq_line()

## Test homoscedasticity (variance assumption)
ggplot(shn_tbl, aes(site, shannon, group = treatment, color = treatment)) +
  geom_jitter()
ggplot(smp_tbl, aes(site, simpson, group = treatment, color = treatment)) +
  geom_jitter()
leveneTest(shannon ~ site, data = shn_tbl) #p = 0.715, null accepted
leveneTest(simpson ~ site, data = smp_tbl) #p = 0.6645, null accepted
plot(shn_tbl$shannon)                      # plot variance, widens
plot(smp_tbl$simpson)                      # plot variance, widens
boxplot(shn_tbl$shannon)                   # plot spread, not symmetrical
boxplot(smp_tbl$simpson)                   # plot spread, not symmetrical
ggplot(shn_tbl, aes(x = shannon,           # plot density, a mess
                    color = site,
                    fill = site)) +
  geom_density(alpha = 0.4)
ggplot(smp_tbl, aes(x = simpson,           # plot density, a mess that's skewed to the right
                    color = site,
                    fill = site)) +
  geom_density(alpha = 0.4)

## NONPARAMETRIC Wilcoxon signed-rank test
#wilcox: Is there a difference between the summed ranks of the differences in ASV counts for the treatments?
shn_ptable <- shn_tbl |>                                      # use count table
  summarise(.by = c(treatment, pair),                         # group data by treatment and pair
            shannon = mean(shannon)) |>                       # count number of unique ASV for each treatment*pair
  pivot_wider(names_from = treatment, values_from = shannon)  # seperate columns for treatment counts
shn_ptable
wilcox.test(shn_ptable$NoRxF, shn_ptable$RxF, paired = TRUE)  #p = 0.5625, null accepted
smp_ptable <- smp_tbl |>                                      # use count table
  summarise(.by = c(treatment, pair),                         # group data by treatment and pair
            simpson = mean(simpson)) |>                       # count number of unique ASV for each treatment*pair
  pivot_wider(names_from = treatment, values_from = simpson)  # seperate columns for treatment counts
smp_ptable
wilcox.test(smp_ptable$NoRxF, smp_ptable$RxF, paired = TRUE)  #p = 0.8438, null accepted

## NONPARAMETRIC Wilcoxon signed-rank test
#wilcox: Is there a difference between the summed ranks of the differences in ASV counts for the months?
shn_ptable2 <- shn_tbl |>                                # use count table
  summarise(.by = c(month, pair),                        # group data by treatment and pair
            shannon = mean(shannon)) |>                  # count number of unique ASV for each treatment*pair
  pivot_wider(names_from = month, values_from = shannon) # seperate columns for treatment counts
shn_ptable2
wilcox.test(shn_ptable2$Jan, shn_ptable2$July, paired = TRUE)  #p = 0.5625, null accepted
smp_ptable2 <- smp_tbl |>                                      # use count table
  summarise(.by = c(month, pair),                              # group data by treatment and pair
            simpson = mean(simpson)) |>                        # count number of unique ASV for each treatment*pair
  pivot_wider(names_from = month, values_from = simpson)       # seperate columns for treatment counts
smp_ptable2
wilcox.test(smp_ptable2$Jan, smp_ptable2$July, paired = TRUE)  #p = 0.4375, null accepted

## NONPARAMETRIC Friedman Test
#friedman: Is there a difference between the treatment rank sums?
shn_gtable <- shn_tbl |>
  summarise(.by = c(treatment, month, pair),    # group data by treatment and month
            avg_shn = mean(shannon)) |>         # count number of unique ASV for each treatment*month
  pivot_wider(names_from = c(treatment, month), # seperate columns for treatment*month counts
              values_from = avg_shn)
shn_gmtrx <- as.matrix(shn_gtable[, 2:5])       # convert to matrix for Friedman Test
shn_gmtrx
friedman.test(shn_gmtrx)                        #p = 0.7055, null accepted
smp_gtable <- smp_tbl |>
  summarise(.by = c(treatment, month, pair),    # group data by treatment and month
            avg_smp = mean(simpson)) |>         # count number of unique ASV for each treatment*month
  pivot_wider(names_from = c(treatment, month), # seperate columns for treatment*month counts
              values_from = avg_smp)
smp_gmtrx <- as.matrix(smp_gtable[, 2:5])       # convert to matrix for Friedman Test
smp_gmtrx
friedman.test(smp_gmtrx)                        #p = 0.9776, null accepted

## NOTE:
# There is no evidence to support that simpson or shannon diversity differ across treatments, timepoints,
# or between treatments across timepoints.


#----------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------- BETA DIVERSITY; NMDS PLOT ----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

## Create freq matrix
frq_mat <- read.delim("data/rarified_feat_freq_full.tsv", sep = "\t", check.names = FALSE, skip = 1) |>
  as.data.frame() |>
  rename_with(~ paste0("sample_", .x, recycle0 = TRUE)) |>      # rename columns to not begin w/ integer
  rename(ASV = "sample_#OTU ID") |>                             # rename ASV column to ASV
  pivot_longer(-ASV) |>
  mutate(value = as.integer(value)) |>
  rename(subsite = "name", name = "ASV") |>                     # rename ASV column to ASV
  group_by(name) |>
  mutate(total = sum(value)) |>
  filter(total != 0) |>
  ungroup() |>
  select(subsite, name, value) |>
  pivot_wider(id_cols = subsite) |>
  column_to_rownames("subsite")
frq_mat <- as.matrix(frq_mat)
frq_mat[1:20, 1:2]

## Generate Bray Curtis distance matrix
bray <- vegdist(frq_mat, method = "bray")
str(bray)                         # look at structure
samples <- attr(bray, "Labels")   # access labels
bray_tbl <- as_tibble(bray) |>
  cbind(samples)
bray_tbl

## Create function to extract metadata from sample name
sample_lookup <- function(df, sampleid) {
  df |>
    mutate(site = factor(str_extract(df[[sampleid]], "(?<=_)\\d+"),    # extract after "_" through digits
                         levels = c("2", "1", "4", "3",
                                    "5", "6", "7", "8",
                                    "10", "9", "12", "11")),           # order sites, NoRxF before RxF
           subsite = str_extract(df[[sampleid]], "(?<=_)\\d+\\w"),     # after "_" through digits and one word char
           month = str_replace_all(str_extract(df[[sampleid]], "(\\d+)$"),     # all digits from end of string
                                   c("1" = "January", "7" = "July")),                 # format months
           treatment = factor(str_replace_all(site, c("10" = "NoRxF", "9" = "RxF",    # designate treatments
                                                      "12" = "NoRxF", "11" = "RxF",   # double-digit #'s first
                                                      "2" = "NoRxF", "1" = "RxF",
                                                      "4" = "NoRxF", "3" = "RxF",
                                                      "5" = "NoRxF", "6" = "RxF",
                                                      "7" = "NoRxF", "8" = "RxF"))),
           pair = str_replace_all(site, c("10" = "PairY", "9" = "PairY",              # designate pairs
                                          "12" = "PairZ", "11" = "PairZ",             # double-digit #'s first
                                          "2" = "PairT", "1" = "PairT",
                                          "4" = "PairV", "3" = "PairV",
                                          "5" = "PairW", "6" = "PairW",
                                          "7" = "PairX", "8" = "PairX")),
           tmntmonth = factor(paste0(treatment, month),
                              levels = c("NoRxFJanuary", "NoRxFJuly",
                                         "RxFJanuary", "RxFJuly"),
                              labels = c("Wildfire Only, January",
                                         "Wildfire Only, July",
                                         "RxF + Wildfire, January",
                                         "RxF + Wildfire, July")))
}

## NMDS ordination using Bray-Curtis distance matrix
set.seed(10192026)
nmds <- metaMDS(bray)
captit <- "Matrix Used: <i>BrayCurtis</i>,Stress: <i>0.2119787</i>"
plttit <- "NMDS Ordination Plot for Bray Curtis Distance"
ccap <- "#363535"

## Plot NDMS ordination colored by treatment
scores(nmds) |>   # by treatment
  as_tibble(rownames = "samples") |>
  sample_lookup(sampleid = "samples") |>
  ggplot(aes(x = NMDS1, y = NMDS2, color = treatment)) +
  geom_point() +
  labs(title = plttit,
       caption = captit) +
  theme(plot.title.position = "plot",                         # align title to plot
        plot.title = element_text(face = "bold"),             # make title bold
        plot.caption = element_markdown(hjust = 0, color = ccap), # move caption to be left-justified
        plot.caption.position = "plot")                       # align caption to plot

## Plot NDMS ordination colored by month
scores(nmds) |>   # by month
  as_tibble(rownames = "samples") |>
  sample_lookup(sampleid = "samples") |>
  ggplot(aes(x = NMDS1, y = NMDS2, color = month)) +
  geom_point() +
  labs(title = plttit,
       caption = captit) +
  theme(plot.title.position = "plot",                         # align title to plot
        plot.title = element_text(face = "bold"),             # make title bold
        plot.caption = element_markdown(hjust = 0, color = ccap), # move caption to be left-justified
        plot.caption.position = "plot")                       # align caption to plot

## Plot NDMS ordination colored by treatment and month
scores(nmds) |>   # by tmntmonth
  as_tibble(rownames = "samples") |>
  sample_lookup(sampleid = "samples") |>
  ggplot(aes(x = NMDS1, y = NMDS2, color = tmntmonth)) +
  geom_point() +
  geom_text(aes(label = subsite), position = "nudge") +    # label point w/ site
  scale_color_manual(values = c("#202020", "#a7a3a3", "#ac2b2b", "#f38383"),
                     name = NULL) +
  labs(title = plttit,
       caption = captit) +
  theme_minimal() +
  theme(plot.title.position = "plot",                             # align title to plot
        plot.title = element_text(face = "bold"),                 # make title bold
        plot.caption = element_markdown(hjust = 0, color = ccap), # move caption to be left-justified
        plot.caption.position = "plot",                           # align caption to plot
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

## Plot NDMS ordination by site
scores(nmds) |>   # by site
  as_tibble(rownames = "samples") |>
  sample_lookup(sampleid = "samples") |>
  ggplot(aes(x = NMDS1, y = NMDS2, color = site)) +
  geom_point() +
  scale_color_brewer(palette = "Paired") +
  labs(title = plttit,
       caption = captit) +
  theme(plot.title.position = "plot",                         # align title to plot
        plot.title = element_text(face = "bold"),             # make title bold
        plot.caption = element_markdown(hjust = 0, color = ccap), # move caption to be left-justified
        plot.caption.position = "plot")                       # align caption to plot


#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------- BETA DIVERSITY; BIPLOT -----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#inspired by https://www.youtube.com/watch?v=_EjNR6YdnGM&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=170

## Get NMDS scores and add metadata
nmds_pos <- scores(nmds) |>           # use scores from NMDS section above
  as_tibble(rownames = "samples") |>  # convert to tibble df
  sample_lookup(sampleid = "samples") # use function defined in NMDS section

## Get the correlation value for abundance of ASV and NMDS x-coordinate
nmds_cor_x <- read.delim("data/rarified_feat_freq_full.tsv", sep = "\t", check.names = FALSE, skip = 1) |>
  as.data.frame() |>
  rename_with(~ paste0("sample_", .x, recycle0 = TRUE)) |>        # rename columns to not begin w/ integer
  rename(asv = "sample_#OTU ID") |>                               # rename ASV column to asv
  pivot_longer(-asv, values_to = "abundance", names_to = "samples") |>
  rename_all(tolower) |>                                          # make all columns lowercase
  inner_join(nmds_pos, by = "samples") |>
  nest(data = -asv) |>
  mutate(cor_x = map(data, ~cor.test(.x$abundance,
                                     .x$NMDS1,
                                     method = "spearman",
                                     exact = FALSE) |>
                       tidy())) |>
  unnest(cor_x) |>
  select(asv, estimate, p.value)

## Get the correlation value for abundance of ASV and NMDS y-coordinate
nmds_cor_y <- read.delim("data/rarified_feat_freq_full.tsv", sep = "\t", check.names = FALSE, skip = 1) |>
  as.data.frame() |>
  rename_with(~ paste0("sample_", .x, recycle0 = TRUE)) |>        # rename columns to not begin w/ integer
  rename(asv = "sample_#OTU ID") |>                               # rename ASV column to asv
  pivot_longer(-asv, values_to = "abundance", names_to = "samples") |>
  rename_all(tolower) |>                                          # make all columns lowercase
  inner_join(nmds_pos, by = "samples") |>
  nest(data = -asv) |>
  mutate(cor_x = map(data, ~cor.test(.x$abundance,
                                     .x$NMDS2,
                                     method = "spearman",
                                     exact = FALSE) |>
                       tidy())) |>
  unnest(cor_x) |>
  select(asv, estimate, p.value)

## Join NMDS correlation x and y values
nmds_cors <- inner_join(nmds_cor_x, nmds_cor_y, by = "asv") |> #2,902 rows
  #filter(p.value.x < 0.05 | p.value.y < 0.05) |> #417 rows
  #filter(p.value.x < 0.01 | p.value.y < 0.01) |> #171 rows
  filter(abs(estimate.x) > 0.5 | abs(estimate.y) > 0.5) #5 rows

## Plot the ASV correlation vectors over ordination
nmds_pos |>
  ggplot(aes(x = NMDS1, y = NMDS2, color = tmntmonth)) +
  geom_point() +
  geom_text(aes(label = subsite), position = "nudge") +    # label point w/ site
  scale_color_manual(values = c("#202020", "#a7a3a3", "#ac2b2b", "#f38383"),
                     name = NULL) +
  labs(title = plttit,
       caption = captit) +
  theme_minimal() +
  theme(plot.title.position = "plot",                             # align title to plot
        plot.title = element_text(face = "bold"),                 # make title bold
        plot.caption = element_markdown(hjust = 0, color = ccap), # move caption to be left-justified
        plot.caption.position = "plot",                           # align caption to plot
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  geom_segment(data = nmds_cors,
               aes(x = 0, xend = estimate.x, y = 0, yend = estimate.y),
               inherit.aes = FALSE) +
  geom_text_repel(data = nmds_cors,
                  aes(x = estimate.x, y = estimate.y, label = asv),
                  inherit.aes = FALSE)


#----------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------- CHANGE RELATIVE TO POINT ----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
#inspired by https://www.youtube.com/watch?v=jLVKJ_n6Qd0&list=PLmNrK_nkqBpKY3SZiivlIGvcLX-KHmfR8&index=167

bray |>                                                # use distance matrix from NMDS section
  as.matrix() |>                                       # convert to matrix
  as_tibble(rownames = "sample") |>                    # convert to tibble
  pivot_longer(-sample) |>                             # make tidy
  filter(sample < name) |>                             # filter for only one triangle of matrix
  sample_lookup(sampleid = "sample") |>                # use function defined in NMDS section
  rename(site_a = "site",
         subsite_a = "subsite",
         month_a = "month",
         treatment_a = "treatment",
         pair_a = "pair") |>
  sample_lookup(sampleid = "name") |>                  # use function defined in NMDS section
  select(-tmntmonth) |>
  #filter(month_a == month & treatment_a != treatment)  # across treatments
  filter(subsite_a == subsite) |>                      # across months
  ggplot(aes(x = subsite, y = value, group = site)) +
  geom_line()
