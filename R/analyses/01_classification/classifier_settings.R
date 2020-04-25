## PREAMBLE ===================================================================
# script what to run, already set like in paper;
# tip: set runs to a smaller number than 1000 if you do not want to wait hours
# number of runs to get the CV results distribution, >=1000 recommended
# runs also will be the name of the results folder
# available results in results folder and their contents described further down
# for a result (try runs = 10 for starters)

## WHAT TO RUN ================================================================
# just the behavioral parameter sets
outer_cv_noaddfeat = F # with outer CV, getting generalization error, Ha
noout_cv_noaddfeat = F # no outer CV, get complete model on whole sample

# only MRI features (but cut for only cue reactivity (CR))
outer_cv_addfeaton = T # Ha only, i.e. MRI only
noout_cv_addfeaton = T # to get the model estimation on whole data set 

# behavior plus MRI
outer_cv_wiaddfeat = F # behavior plus physio, Ha (MRI)
noout_cv_wiaddfeat = F # behavior plus physio, get complete model

# control model
outer_cv_c_model = F # baseline model for classification; predict with covariate only

# cut fmri data for cue reactivity (cr) only?
cut_for_mri_cr_only = F # the fmri data will be cut for cr data only

# report
do_report_feat_only = T
do_report_with_added_feat = F

# number of runs to get the CV results distribution, >=1000 recommended
# runs also will be the name of the results folder
# 1000: fMRI against control model (control is baseline model, i.e. some control variable)
# 350: fMRI with behavioral features
# 351: only behavioral features
# 352: fMRI only CR
# 353: behavioral plus fmri (only CR)
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

