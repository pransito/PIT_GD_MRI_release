

# report on performance for multiple classification attempts
if (init_done == F) {
  stop('You have to first run the select_study.R script with which_study set to "MRI" ')
}

stopifnot(which_study == 'MRI') 

# get the functions
source('report_mult_classif_funcs.R')

# get all the classifiers and their performances in AUC
res = get_aucs_of_classifiers()

# get the random classification
random_class_unb = get_random_classification(runs=4000, direction='<')
random_class_b = get_random_classification(runs=4000, direction='auto')

# get the performance of old classifier
aucs_old_class = get_performance_of_prior_study_classif(all_aucs=random_class_unb$all_aucs, 
                                                        all_accs=random_class_unb$all_accs,
                                                        all_sens=random_class_unb$all_sens,
                                                        all_spec=random_class_unb$all_spec)


# bind the random classification on top and the old classifier to the bottom
cur_random_unb = list(random_unbiased = random_class_unb$all_aucs)
cur_random_b = list(random_biased = random_class_b$all_aucs)
#cur_old_class = list(behav_Genauck_2019 = aucs_old_class$real_aucs)
res = c(cur_random_unb, cur_random_b, res)

# make a df for density plotting
df_performance_classifiers = make_df_performance_classifiers(res)
names(df_performance_classifiers) = c('classifier', 'auc')

# get the data-frame of performance
perf = aggregate(df_performance_classifiers$auc, by=list(df_performance_classifiers$classifier), 
                 FUN=agk.mean.quantile.c, lower=0.095, upper=0.975)
perf_gr=perf$Group.1
perf=data.frame(perf$x)
perf$classifier = perf_gr
perf = perf[order(perf$mean),]
perf$classifier = factor(as.character(perf$classifier), levels=perf$classifier)
df_performance_classifiers$classifier = factor(df_performance_classifiers$classifier, 
                                               levels = levels(perf$classifier))
perf[c('mean', 'lower', 'upper')] = round(perf[c('mean', 'lower', 'upper')],2)
write.table(perf, file='results/classification_performance.csv')

# plot
p = make_density_plots(cur_dat = df_performance_classifiers)

# p value fmri only against fmri+behav
fmri_with_behav = df_performance_classifiers$auc[df_performance_classifiers$classifier=='fmri_with_behav']
behav_only = df_performance_classifiers$auc[df_performance_classifiers$classifier=='behav_features_only']
agk.density_p.c(fmri_with_behav-behav_only, 0)

