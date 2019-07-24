R Code for the PIT Gambling Disorder paper (MRI data)
-----------------------------------------------------------

Code for publication of "Neural correlates of cue-induced changes in decision-making distinguish subjects with gambling disorder from healthy controls"
by: Alexander Genauck (2019)


How to use
----------

Fork the whole repository or download the zip and put it into some working directory.
The .RData file has all the data and functions you will need. Especially data_pdt (choice data) and dat_match (questionnaire, demographic data). You will need R (https://cran.r-project.org/bin/windows/base/) and Rstudio (https://www.rstudio.com/products/rstudio/#Desktop) for this. Both are freely available software.

1)
Run the "R/select_study.R" script. It selects the Cohort (MRI) and initializes all the analyses. Do this before anything else. 

Careful: The script "R/analyses/01_classification/group_pred_loop_v7.R" (called by "R/select_study.R") installs and loads many R packages (see in the beginning of that script). They are all useful and should not hurt your installation of R or bother your other code. However, revise the list before you run the code and decide
if you would like to continue.

2)
If you want to see results you set everything to FALSE in "WHAT TO RUN" section in "R/analyses/01_classification/classifier_settings.R", except "do_report_feat_only" to TRUE. Set runs to "1000" to see the results reported in the paper, then run the script: R/analyses/01_classification/reporting_v7.R

3)
The machine learning part is started with the script "R/analyses/01_classification/group_pred_loop_v7.R" Before running it, you must adjust the settings in: "R/analyses/01_classification/classifier_settings.R" [set "runs" to a low number like 10 for it is an intense script which takes a long time to run].
Set outer_cv_addfeaton = 1 and noout_cv_addfeaton = 1 and do_report_feat_only = 0 (unless "runs" is > 30)

4)
The classical group analysis part using hierarchical regression (lme4) is  with the
"/02_classical_group_analyses/glmer_accRate_la_cat_v3.R script"; Before running it, make sure to run "R/select_study.R" first with which_study set to "MRI". It is set such that the glmer models are run (10 to 20 minutes) and compared according to main text and supplements
a plot is generated showing the shifts in acceptance rate according to cue category

5)
Making the ratings graph: is automatically done when running "R/select_study.R" (see)

6)
Ratings: statistical tests. Run the script "R/analyses/03_image_adequacy/ratings_analysis_for_paper.R". Check the instructions at the top of the script.

7)
The Matlab folder holds information on the sequences used for functional magnetic resonance imaging and also the SPM12 batch for preprocessing and the single subject model.


Enjoy! Any questions, or bugs please send by raising an issue in this GitHub repository.