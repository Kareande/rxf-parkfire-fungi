# KAT ANDERSON AVITI SEQUENCE ANALYSIS WITH R
# Last updated 2026: 

#-------------------------------------------------------------#
#-------------------- PREPPING WORKSPACE ---------------------#
#-------------------------------------------------------------#

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
library(rstatix)        # stats package
library(ggpubr)         # idk what this was for tbh

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

#------------------------------------------#
#---- GENERATE COLOR PALLET FOR PLOTS -----#
#------------------------------------------#

## Set seed for colors (repeatability)
seed <- 82526
set.seed(seed)

## Sample n colors from ColorBrewer
n <- 74
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                            rownames(qual_col_pals)))

pie(rep(1, n), col = col_vector, main = paste0("Color Pallet ", seed, "
(',n,' Colors)"))

col_vector <- sort(col_vector)
n_colors <- sample(col_vector, n)

## Generate pallet of n colors
par(mfrow = c(1, 1), bg = "white")
pie(rep(1, n), col = n_colors, main = paste0("Color Pallet ", seed, "
(',n,' Colors)"))

## Save generated pallet
dev.print(device = png,
          filename = paste0("Images/ColorPallet_", seed, "_", n, "Col.png"),
          width = 672, height = 672)
dev.off()


#-------------------------------------------------------------#
#--------------------- TIDY DATASET --------------------------#
#-------------------------------------------------------------#
## NOTE:
# I filtered out all samples with less than 500 reads,
## Load in data
frq_tbl <- read.delim("rarified_feat_freq_full.tsv", sep = "\t", check.names = FALSE, skip = 1) |>
  as.data.frame() |>
  rename_with(~ paste0("sample_", .x, recycle0 = TRUE)) |>      # rename columns to not begin w/ integer
  rename(ASV = "sample_#OTU ID") |>                             # rename ASV column to ASV
  pivot_longer(-ASV, values_to = "Count") |>                    # make tidy with one row per ASV-Sample
  mutate(Site = substr(name, 8, nchar(name) - 2),                 # extract Site from name (ex: 10 from sample_10B1)
         Subsite = substr(name, 8, nchar(name) - 1)) |>           # extract Subsite from name (ex: 10B sample_10B1)
  mutate(CollectMonth = substr(name, nchar(name), nchar(name))) |>                         # extract collection month
  mutate(CollectMonth = str_replace_all(CollectMonth, c("1" = "Jan", "7" = "July"))) |>    # reformat months
  mutate(Treatment = str_replace_all(Site, c("10" = "NoRxF", "11" = "RxF", "12" = "NoRxF", # reformat treatments
                                             "1" = "RxF", "2" = "NoRxF", "3" = "RxF",
                                             "4" = "NoRxF", "5" = "NoRxF", "6" = "RxF",
                                             "7" = "NoRxF", "8" = "RxF", "9" = "RxF"))) |>
  mutate(Pair = str_replace_all(Site, c("10" = "PairY", "11" = "PairZ", "12" = "PairZ",    # designate pairs
                                        "1" = "PairT", "2" = "PairT", "3" = "PairV",
                                        "4" = "PairV", "5" = "PairW", "6" = "PairW",
                                        "7" = "PairX", "8" = "PairX", "9" = "PairY"))) |>
  mutate(TmntPair = factor(Site,                                                           # create factors
                           levels = c("2", "1", "4", "3", "5", "6", "7", "8", "10", "9", "12", "11"),
                           labels = c("NoRxF Pair T", "RxF Pair T", "NoRxF Pair V", "RxF Pair V",
                                      "NoRxF Pair W", "RxF Pair W", "NoRxF Pair X", "RxF Pair X",
                                      "NoRxF Pair Y", "RxF Pair Y", "NoRxF Pair Z", "RxF Pair Z"))) |>
  mutate(Treatment = factor(Treatment),
         Pair = factor(Pair),
         Site = factor(Site, levels = c("2", "1", "4", "3", "5", "6",                      # RxF before NoRxF site
                                        "7", "8", "10", "9", "12", "11"))) |>
  select(ASV, Site, Subsite, Treatment, CollectMonth, Pair, TmntPair, Count)               # re-order desired columns
frq_tbl

## Histogram
ggplot(frq_tbl, aes(x = Count)) +
  geom_histogram(fill = n_colors[33], color = n_colors[35], binwidth = 20)  # yikes

## Remove any ASVs with zero counts
frq_tbl <- filter(frq_tbl, Count > 0)

## Histogram
ggplot(frq_tbl, aes(x = Count)) +
  geom_histogram(fill = n_colors[33], color = n_colors[35], binwidth = 20)  # somewhat better

## Visualize median observations per treatment pair over across timepoints
frq_tbl |> mutate(CollectMonth = factor(CollectMonth,
                                        levels = c("Jan", "July"),
                                        labels = c("Samples Collected 6mo Post-fire",
                                                  "Samples Collected 12mo Post-fire"))) |>
  ggplot(aes(x = Pair, y = Count, color = Treatment)) +
  facet_wrap(~ CollectMonth, nrow = 1, scales = "free_x", strip.position = "bottom") +
  stat_summary(fun.data = median_hilow,
               fun.args = list(conf.int = 0.5),
               geom = "pointrange",
               position = position_dodge(width = 0.5)) +
  theme(strip.placement = "outside",
        strip.background = element_blank()) +
  labs(x = NULL,
       y = "ASV Count") +
  scale_color_manual(name = NULL,
                     breaks = c("NoRxF", "RxF"),
                     labels = c("Wildfire Only", "RxF + Wildfire"),
                     values = c("black", "red"))

## Mean observations per subsite/sample (total abundance is 6222 sampling depth)
avg_tbl <- frq_tbl |>
  summarise(.by = c(Subsite, Treatment, CollectMonth, Pair, TmntPair), # group data
            Avg_ASV = mean(Count))                                           # avg observations
avg_tbl

## Mean observations per site
avg_tbl_site <- frq_tbl |>
  summarise(.by = c(Site, Treatment, CollectMonth, Pair, TmntPair),  # group data
            Avg_ASV = mean(Count))                                   # avg observations
avg_tbl_site

## Visualize mean ASV counts per paired site
asv_plot <- avg_tbl |>
  summarise(.by = c(Pair, Treatment), Avg_ASV = mean(Avg_ASV)) |>
  ggplot(aes(x = Pair, y = Avg_ASV, fill = Treatment)) +
  geom_col(position = "dodge", alpha = 0.75) +
  scale_fill_manual(values = c(col_vector[16], col_vector[22])) +
  labs(title = "Mean ASV Abundance per Paired Site",
       x = "Site Pair",
       y = "Mean ASV Abundance",
       fill = "Treatment") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 20))
asv_plot
ggsave("C://Users//Kat//VS_Code//rxf-parkfire-fungi//figures//Avg_ASV_Pair.png")

## Visualize mean ASV counts per paired site across timepoints
pasv_plot <- avg_tbl_site |>
  mutate(CollectMonth = factor(CollectMonth,
                               levels = c("Jan", "July"),
                               labels = c("Samples Collected 6mo Post-fire",
                                          "Samples Collected 12mo Post-fire"))) |>
  ggplot(aes(x = Pair, y = Avg_ASV, fill = Treatment)) +
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
ggsave("C://Users//Kat//VS_Code//rxf-parkfire-fungi//figures//Avg_ASV_Pair_Mo.png")


#-------------------------------------------------------------------------------------------#
#---------- DOES ASV RICHNESS DIFFER B/T TREATMENTS, MONTHS, OR TREATMENT*MONTHS? ----------#
#-------------------------------------------------------------------------------------------#
#H0: The AVS counts do not differ across treatments
#H1: The AVS counts are different for soil treated and not treated with RxF

## Count unique ASVs per subsite/sample
rich_tbl <- frq_tbl |>                                                      # use freq table
  summarise(.by = c(Site, Subsite, Treatment,                               # group data
                    Pair, TmntPair, CollectMonth),
            count_asv = length(unique(ASV)))                                # number of unique ASV per subsite
rich_tbl

## Check for outliers (keep running until all identified/removed)
identify_outliers(rich_tbl, count_asv)
rich_tbl <- filter(rich_tbl, Subsite != "3E")  #had to remove 3E, 12R, 3D
rich_tbl <- filter(rich_tbl, Subsite != "12R") #had to remove 3E, 12R, 3D
rich_tbl <- filter(rich_tbl, Subsite != "3D")  #had to remove 3E, 12R, 3D

## Test normality (distribution assumption)
shapiro.test(rich_tbl$count_asv)            #p = 0.5808, null accepted
ks.test(rich_tbl$count_asv, "pnorm")        #p < 2.2e-16, null rejected
ad.test(rich_tbl$count_asv)                 #p = 0.4196, null accepted
plot(density(rich_tbl$count_asv))           # plot density, normalish
ggplot(rich_tbl, aes(sample = count_asv)) + # plot Q-Q plot, almost linear
  geom_qq() +
  geom_qq_line()

## Test homoscedasticity (variance assumption)
ggplot(rich_tbl, aes(Site, count_asv, group = Treatment, color = Treatment)) +
  geom_jitter()
leveneTest(count_asv ~ Site, data = rich_tbl) #p = 0.7731, null accepted
plot(rich_tbl$count_asv)                      # plot variance, fine
boxplot(rich_tbl$count_asv)                   # plot spread, almost symmetrical

## PARAMETRIC Paired t-Test, NONPARAMETRIC Wilcoxon signed-rank test
#t.test: Is there is a difference between the mean differences in ASV counts for the Treatments?
#wilcox: Is there a difference between the summed ranks of the differences in ASV counts for the Treatments?
count_ptable <- rich_tbl |>                                      # use count table
  summarise(.by = c(Treatment, Pair),                            # group data by Tmnt and Pair
            count_asv = sum(count_asv)) |>                       # count number of unique ASV for each Tmnt*Pair
  pivot_wider(names_from = Treatment, values_from = count_asv)   # seperate columns for Tmnt counts
count_ptable
t.test(count_ptable$NoRxF, count_ptable$RxF, paired = TRUE)      #p = 0.8132, null accepted
wilcox.test(count_ptable$NoRxF, count_ptable$RxF, paired = TRUE) #p = 0.6875, null accepted

## PARAMETRIC Paired t-Test, NONPARAMETRIC Wilcoxon signed-rank test
#t.test: Is there is a difference between the mean differences in ASV counts for the Months?
#wilcox: Is there a difference between the summed ranks of the differences in ASV counts for the Months?
count_ptable2 <- rich_tbl |>                                      # use count table
  summarise(.by = c(CollectMonth, Pair),                          # group data by Tmnt and Pair
            count_asv = sum(count_asv)) |>                        # count number of unique ASV for each Tmnt*Pair
  pivot_wider(names_from = CollectMonth, values_from = count_asv) # seperate columns for Tmnt counts
count_ptable2
t.test(count_ptable2$Jan, count_ptable2$July, paired = TRUE)      #p = 0.5118, null accepted
wilcox.test(count_ptable2$Jan, count_ptable2$July, paired = TRUE) #p = 0.5625, null accepted

## PARAMETRIC ANOVA repeated measures, NONPARAMETRIC Friedman Test
#no benefit in choosing t.test vs ANOVA
#anova: Is there a difference between the mean ASV counts for the treatments?
#friedman: Is there a difference between the treatment rank sums?
count_gtable <- rich_tbl |>                         # use count table
  summarise(.by = c(Treatment, CollectMonth, Pair), # group data by Treatment and Month
            count_asv = mean(count_asv))            # count number of unique ASV for each Tmnt*Month
print(count_gtable, n = 24)
aov_res <- anova_test(data = count_gtable,          # fit a model, ASV richness in response to Tmnt
                      dv = count_asv,               # dependant variable is count_asv
                      wid = Pair,                   # subject identified is Site
                      within = c(CollectMonth, Treatment))
get_anova_table(aov_res)                            #p > 0.05, null accepted
count_gtbl <- count_gtable |>                           # convert to matrix for Friedman Test
  pivot_wider(names_from = c(Treatment, CollectMonth),   # seperate columns for Tmnt*Month counts
              values_from = count_asv)
count_gmtrx <- as.matrix(count_gtbl[, 2:5])             # convert to matrix for Friedman Test
count_gmtrx
friedman.test(count_gmtrx)                               #p = 0.7055, null accepted

## NOTE:
# Data is normal enough to justify using parametric tests. Did both, but will only use parametric results.
# There is no support for the hypothesis that ASV richness is different between treatments, timepoints, nor
# treatments across timepoints.

## Test timepoint difference Roualdes' way
## Does the mean differences in ASV counts significantly differ across treatments (all sites) or pairs?
diff_table <- count_gtbl |>                              # use Friedman Test table
  mutate(NoRxF_diff = NoRxF_July - NoRxF_Jan, Pair,    # subtract Jan counts from July
         RxF_diff = RxF_July - RxF_Jan) |>
  select(Pair, RxF_diff, NoRxF_diff) |>
  pivot_longer(c(RxF_diff, NoRxF_diff),
               names_to = "Tmnt",
               values_to = "diff_asv") |>
  mutate(Tmnt = substr(Tmnt, 1, nchar(Tmnt) - 5)) # trim the site from subsite
diff_table

## PARAMETRIC ANOVA, NONPARAMETRIC Kruskal-Wallis
#anova test sig diff between the means by analyzing variance
fit <- lm(diff_asv ~ Tmnt, data = diff_table)    # fit model using difference table, difference in response to tretament
anova(fit)                                       #p = 0.03106, null rejected

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
# Professor Roualdes initially grouped by Site (summing the counts for each subsite),
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


#-------------------------------------------------------------------------------------------#
#------- DOES SHANNON/SIMPSON DIFFER B/T TREATMENTS, MONTHS, OR TREATMENT*MONTHS? ----------#
#-------------------------------------------------------------------------------------------#

## Calculate diversity indexes
shn_tbl <- frq_tbl |>                                                  # use ASV table
  summarise(.by = c(Site, Pair, Subsite, Treatment, CollectMonth),             # group data
            shannon = diversity(Count, index = "shannon"))
shn_tbl
smp_tbl <- frq_tbl |>                                                  # use ASV table
  summarise(.by = c(Site, Pair, Subsite, Treatment, CollectMonth),             # group data
            simpson = diversity(Count, index = "simpson"))
smp_tbl

## Check for outliers (keep running until all identified/removed)
identify_outliers(shn_tbl, shannon)
identify_outliers(smp_tbl, simpson)
smp_tbl <- filter(smp_tbl, Subsite != "6R") #had to remove 6R, 5G
smp_tbl <- filter(smp_tbl, Subsite != "5G") #had to remove 6R, 5G

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
ggplot(shn_tbl, aes(Site, shannon, group = Treatment, color = Treatment)) +
  geom_jitter()
ggplot(smp_tbl, aes(Site, simpson, group = Treatment, color = Treatment)) +
  geom_jitter()
leveneTest(shannon ~ Site, data = shn_tbl) #p = 0.715, null accepted
leveneTest(simpson ~ Site, data = smp_tbl) #p = 0.6645, null accepted
plot(shn_tbl$shannon)                      # plot variance, widens
plot(smp_tbl$simpson)                      # plot variance, widens
boxplot(shn_tbl$shannon)                   # plot spread, not symmetrical
boxplot(smp_tbl$simpson)                   # plot spread, not symmetrical
ggplot(shn_tbl, aes(x = shannon,           # plot density, a mess
                    color = Site,
                    fill = Site)) +
  geom_density(alpha = 0.4)
ggplot(smp_tbl, aes(x = simpson,           # plot density, a mess that's skewed to the right
                    color = Site,
                    fill = Site)) +
  geom_density(alpha = 0.4)

## NONPARAMETRIC Wilcoxon signed-rank test
#wilcox: Is there a difference between the summed ranks of the differences in ASV counts for the treatments?
shn_ptable <- shn_tbl |>                                      # use count table
  summarise(.by = c(Treatment, Pair),                         # group data by Tmnt and Pair
            shannon = mean(shannon)) |>                       # count number of unique ASV for each Tmnt*Pair
  pivot_wider(names_from = Treatment, values_from = shannon)  # seperate columns for Tmnt counts
shn_ptable
wilcox.test(shn_ptable$NoRxF, shn_ptable$RxF, paired = TRUE)  #p = 0.5625, null accepted
smp_ptable <- smp_tbl |>                                      # use count table
  summarise(.by = c(Treatment, Pair),                         # group data by Tmnt and Pair
            simpson = mean(simpson)) |>                       # count number of unique ASV for each Tmnt*Pair
  pivot_wider(names_from = Treatment, values_from = simpson)  # seperate columns for Tmnt counts
smp_ptable
wilcox.test(smp_ptable$NoRxF, smp_ptable$RxF, paired = TRUE)  #p = 0.8438, null accepted

## NONPARAMETRIC Wilcoxon signed-rank test
#wilcox: Is there a difference between the summed ranks of the differences in ASV counts for the Months?
shn_ptable2 <- shn_tbl |>                                       # use count table
  summarise(.by = c(CollectMonth, Pair),                        # group data by Tmnt and Pair
            shannon = mean(shannon)) |>                         # count number of unique ASV for each Tmnt*Pair
  pivot_wider(names_from = CollectMonth, values_from = shannon) # seperate columns for Tmnt counts
shn_ptable2
wilcox.test(shn_ptable2$Jan, shn_ptable2$July, paired = TRUE)   #p = 0.5625, null accepted
smp_ptable2 <- smp_tbl |>                                       # use count table
  summarise(.by = c(CollectMonth, Pair),                        # group data by Tmnt and Pair
            simpson = mean(simpson)) |>                         # count number of unique ASV for each Tmnt*Pair
  pivot_wider(names_from = CollectMonth, values_from = simpson) # seperate columns for Tmnt counts
smp_ptable2
wilcox.test(smp_ptable2$Jan, smp_ptable2$July, paired = TRUE)  #p = 0.4375, null accepted

## NONPARAMETRIC Friedman Test
#friedman: Is there a difference between the treatment rank sums?
shn_gtable <- shn_tbl |>
  summarise(.by = c(Treatment, CollectMonth, Pair),    # group data by Treatment and Month
            avg_shn = mean(shannon)) |>                # count number of unique ASV for each Tmnt*Month
  pivot_wider(names_from = c(Treatment, CollectMonth), # seperate columns for Tmnt*Month counts
              values_from = avg_shn)
shn_gmtrx <- as.matrix(shn_gtable[, 2:5])              # convert to matrix for Friedman Test
shn_gmtrx
friedman.test(shn_gmtrx)                               #p = 0.7055, null accepted
smp_gtable <- smp_tbl |>
  summarise(.by = c(Treatment, CollectMonth, Pair),    # group data by Treatment and Month
            avg_smp = mean(simpson)) |>                # count number of unique ASV for each Tmnt*Month
  pivot_wider(names_from = c(Treatment, CollectMonth), # seperate columns for Tmnt*Month counts
              values_from = avg_smp)
smp_gmtrx <- as.matrix(smp_gtable[, 2:5])              # convert to matrix for Friedman Test
smp_gmtrx
friedman.test(smp_gmtrx)                               #p = 0.9776, null accepted

## NOTE:
# There is no evidence to support that simpson or shannon diversity differ across treatments, timepoints,
# or between treatments across timepoints.


#----------------------------------------------------------------------------#
#------- BETA DIVERSITY ----------#
#----------------------------------------------------------------------------#

## Create freq matrix
frq_mat <- read.delim("rarified_feat_freq_full.tsv", sep = "\t", check.names = FALSE, skip = 1) |>
  as.data.frame() |>
  rename_with(~ paste0("sample_", .x, recycle0 = TRUE)) |>      # rename columns to not begin w/ integer
  rename(ASV = "sample_#OTU ID") |>                             # rename ASV column to ASV
  pivot_longer(-ASV) |>
  mutate(value = as.integer(value)) |>
  rename(subsite = "name", name = "ASV") |>                             # rename ASV column to ASV
  group_by(name) |>
  mutate(total = sum(value)) |>
  filter(total != 0) |>
  ungroup() |>
  select(subsite, name, value) |>
  pivot_wider(id_cols = subsite) |>
  as.data.frame()
rownames(frq_mat) <- frq_mat$subsite
frq_mat <- frq_mat[, -1]
frq_mat <- as.matrix(frq_mat)
frq_mat[1:20, 1:2]

## Generate distances
set.seed(10192026)
bray <- vegdist(frq_mat, method = "bray")
nmds <- metaMDS(bray)
scores(nmds) |>
  as_tibble(rownames = "subsite") |>
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point()

nmds <- metaMDS(frq_mat, autotransform = FALSE)
scores(nmds) |>
  as_tibble(rownames = "subsite") |>
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point()



