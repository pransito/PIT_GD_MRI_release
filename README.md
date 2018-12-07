R Code for the PIT Gambling Disorder paper (MRI data)
-----------------------------------------------------------

Code for publication of "Neural correlates of cue-induced changes in decision-making distinguish subjects with gambling disorder from healthy controls"
by: Alexander Genauck (2018)


How to use
----------

Fork the whole repository or download the zip and put it into some working directory.
The .RData file has all the data and functions you will need. Especially data_pdt (choice data) and dat_match (questionnaire, demographic data). You will need R (https://cran.r-project.org/bin/windows/base/) and Rstudio (https://www.rstudio.com/products/rstudio/#Desktop) for this. Both are freely available software.

1)
Run the "R/select_study.R" script. It selects the Cohort (MRI) and initializes all the analyses. Do this before anything else. 

Careful: The script "R/analyses/01_classification/group_pred_loop_v7.R" (called by "R/select_study.R") installs and loads many R packages (see in the beginning of the script). They are all useful and should not hurt your installation of R or bother your other code. However, revise the list before you run the code and decide
if you would like to continue.

2)
The machine learning part is started with the script "R/analyses/01_classification/group_pred_loop_v7.R" Before running it, you may adjust the settings in: "R/analyses/01_classification/classifier_settings.R" [set "runs" to a low number like 10 for it is an intense script which takes a long time to run]  If you want to see results you set everything to FALSE, except the reporting section to TRUE (all of them) Set runs to "1000" to see the results reported in the paper. 

3)
The classical group analysis part using hierarchical regression (lme4) is  with the
"/02_classical_group_analyses/glmer_accRate_la_cat_v3.R script"; Before running it, make sure to run "R/select_study.R" first with which_study set to "MRI". It is set such that the glmer models are run (10 to 20 minutes) and compared according to main text and supplements
a plot is generated showing the shifts in acceptance rate according to cue category

4)
Making the ratisoduced when running "group_pred_loop_v7.R". 

5)
Ratings: statistical tests. Run the script "R/analyses/03_image_adequacy/ratings_analysis_for_paper.R". Check the instructions at the top of the script.


Enjoy! Any questions, or bugs please send by raising an issue in this GitHub repository.