## PREAMBLE ===================================================================
# script what to run, already set like in paper;
# tip: set runs to a smaller number than 1010 if you do not want to wait hours
# number of runs to get the CV results distribution, >=1000 recommended
# runs also will be the name of the results folder
# available results in results folder and their contents described further down
# for a result (try 10 for starters)

## WHAT TO RUN ================================================================
# MRI (saved under prefix "phys" in results)
outer_cv_addfeaton      = F # Ha only, i.e. physio/MRI  (with cross-validation)
noout_cv_addfeaton      = F # to get the complete model (no cross-validation)

# control model
outer_cv_c_model        = F # control model/null-model for classification; predict with covariate only

# what to report
do_report               = T
do_report_feat_only     = T

# number of runs to get the CV results distribution, >=1000 recommended
# runs also will be the name of the results folder
# 1000: fMRI predictors
runs = 1000

# advanced settings for other studies =========================================
# leave this section as is
if (which_study == 'MRI') {
  # Any reporting of p-values against null? Set to F if you do that in a separate script.
  report_CV_p = T
} else {
  # Any reporting of p-values against null? Set to F if you do that in a separate script.
  report_CV_p = T
}

# no other features, only behavior
# master add MRI data
add_cr_pp_ma         = T
# master ratings
add_cr_ra_ma         = F

