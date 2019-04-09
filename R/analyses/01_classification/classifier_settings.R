## PREAMBLE ===================================================================
# script what to run, already set like in paper;
# tip: set runs to a smaller number than 1000 if you do not want to wait hours
# number of runs to get the CV results distribution, >=1000 recommended
# runs also will be the name of the results folder
# available results in results folder and their contents described further down
# for a result (try runs = 10 for starters)

## WHAT TO RUN ================================================================
# only MRI features
outer_cv_addfeaton = 0 # Ha only, i.e. MRI
noout_cv_addfeaton = 0 # to get the model estimation on whole data set 

# control model
outer_cv_c_model = 0 # baseline model for classification; predict with covariate only

# report
do_report_feat_only = 1

# number of runs to get the CV results distribution, >=1000 recommended
# runs also will be the name of the results folder
# 1000: fMRI against control model (control is baseline model, i.e. some control variable)
runs = 1000

# advanced settings for (other studies)leave as is) ===========================
# Any reporting of p-values against null? Set to F if you do that in a separate script.
report_CV_p = T

# no other features, only behavior
# master add cue reactivity: peripheral physiology or MRI
if (outer_cv_noaddfeat == T | noout_cv_noaddfeat == T | do_report_no_added_feat == T) {
  add_cr_pp_ma = F
} else {
  add_cr_pp_ma = T
}

# master add cue reactivity: ratings
# should never be done, cause ratings are post-experiment
add_cr_ra_ma = F

