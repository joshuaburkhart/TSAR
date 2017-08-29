# ---
# title: Subchallenge 1 Analysis
# author: Joshua Burkhart
# date: Aug 16, 2017
# ---

# download data
source("../data/Download.R")

# calculate genesets
## all predictor lists from teams that separated predictors by virus
## including teams who only participated in subchallenge 1
source("./per_study/PerStudyAnalysis.R")

## all predictor lists from teams that had both timepoints from at least
## one of the challenges with significant results (sub2 and sub3) with
## less than all the predictors and valid, parsable files
source("./per_timepoint/PerTimepointAnalysis.R")

## all predictors lists fom teams that had both timepoints from both
## of the challenges with significant results (sub2 and sub3) with
## less than all the predictors and valid, parsable files
source("./per_subchallenge/PerSubchallengeAnalysis.R")