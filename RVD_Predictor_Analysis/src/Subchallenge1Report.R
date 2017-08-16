# ---
# title: Subchallenge 1 Analysis
# author: Joshua Burkhart
# date: Aug 16, 2017
# ---

# download data
source("../data/Download.R")

# calculate genesets
source("./per_study/PerStudyGeneSets.R")
source("./per_subchallenge/PerSubchallengeGeneSets.R")
source("./per_timepoint/PerTimepointGeneSets.R")

# generate heatmaps
source("./per_study/PerStudyHeatmap.R")
source("./per_subchallenge/PerSubchallengeHeatmap.R")
source("./per_timepoint/PerTimepointHeatmap.R")

# perform pathway analysis
source("./per_study/PerStudyPA.R")
source("./per_subchallenge/PerSubchallengePA.R")
source("./per_timepoint/PerTimepointPA.R")