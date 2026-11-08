# rxf-parkfire-fungi
This repo explores diversity and community function using fungal amplicon sequence variants (ASV) to test the hypothesis that soil burned in the 2024 Park Fire only has lower soil fungal diversity than soil treated with a single-entry prescribed fire prior to the wildfire.

Soil was collected from Chico State's Big Chico Creek Ecological Reserve, extracted, amplified, and sequenced data by UC Davis’ Technology Core Facility. Sequence data was provided demultiplexed and in FASTQ format. Sequence processing was done using [Qiime2](https://qiime2.org/), and given the size of sequence data, I used the [National Research Platform](https://nrp.ai/) to deploy storage and Kubernetes pod. The following diversity and community analysis was done using R on a local machine.


## File descriptions
**Anderson-Qiime2-SequenceAnalysis.txt** contains my Qiime2 pipeline with the following steps:
- Importing data into Qiime
- Extraction of ITS Region
- Quality filtering and denoising
- Decontaminating using controls
- Training classifier on UNITE data
- Assigning taxonomy
- Exporting data from Qiime

**Anderson-Qiime2-DiversityAnalysis.txt** contains experimental diversity analysis using Qiime2 functionality.

**Diversity-Analysis.R** is my R script for exploring, plotting, and statistically testing the following alpha diversity metrics on ASV, taxa, and functional groups:
- Relative abundance
- Richness
- Shannon and Simpson Diversity Indexes

**Community-Analysis.R** is my R script for exploring, plotting, and statistically testing the the following beta diversity metrics:
- Bray-Curtis and Jaccard Distance
- NMDS
- Differential abundance

Processed sequence data is available in /data, and plots are available in /figures.
