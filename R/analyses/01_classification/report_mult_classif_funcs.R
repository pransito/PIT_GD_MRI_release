


# functions for multiple classificatin attempts ===============================
get_random_classification = function(runs=4000, direction='auto') {
  
  # TODO: the direction option should be an input argument
  #       i.e. biased vs. perfect randomization
  
  message('Getting random classifier performance...')
  message(paste0('...with direction set to ', direction))
  # get the random classification
  
  # functions
  get.truth = function() {
    sample(c(rep('HC',3),rep('PG',3)))
  }
  
  get.truth.4 = function() {
    sample(c(rep('HC',4),rep('PG',4)))
  }
  
  get.truth.2 = function() {
    sample(c(rep('HC',2),rep('PG',2)))
  }
  
  get.truth.1 = function() {
    sample(c('HC','PG'),size = 1)
  }
  
  # set runs of random classification
  runs0 = runs
  
  # under 0
  # pooled
  all_aucs  = c()
  all_aucsl = list()
  all_accs  = c()
  all_accsl = list()
  all_sens  = c()
  all_spec  = c()
  for (ii in 1:runs0) {
    inner_truths = c()
    inner_resps  = c()
    # get 30 in each group
    for (jj in 1:10) {
      # get truth
      inner_truths = c(inner_truths,as.character(get.truth()))
      # get response
      inner_resps  = c(inner_resps,as.numeric(randn(1,6)*10))
    }
    
    # # add 4
    # inner_truths = c(inner_truths,as.character(get.truth.2()))
    # # get response
    # inner_resps  = c(inner_resps,as.numeric(randn(1,4)*10))
    # 
    # # add 1
    # inner_truths = c(inner_truths,as.character(get.truth.1()))
    # # get response
    # inner_resps  = c(inner_resps,as.numeric(randn(1,1)*10))
    
    # cur_auc
    cur_roc         = roc(inner_truths,inner_resps, levels=c('HC', 'PG'), direction=direction)
    all_aucs[ii]    = cur_roc$auc
    all_aucsl[[ii]] = cur_roc
    
    # accuracy
    inner_preds     = ifelse(inner_resps<0,'HC','PG')
    all_accsl[[ii]] = inner_truths == inner_preds
    all_accs[ii]    = mean(all_accsl[[ii]])
    
    # sens and spec
    cur_cm = caret::confusionMatrix(table(inner_truths,inner_preds))
    all_sens[ii]    = cur_cm$byClass[1]
    all_spec[ii]    = cur_cm$byClass[2]
  }
  
  return(list(all_aucs=all_aucs, all_aucsl=all_aucsl, all_accsl=all_accsl, all_accs=all_accs, all_sens=all_sens, 
              all_spec=all_spec))
}


get_performance_of_prior_study_classif = function(all_aucs, all_accs, all_sens, all_spec) {
  
  # function to get ensemble an mean performance of behavioral study (Genauck et al., 2019) classifier on
  # MRI behavioral data
  # the input variables are the performance under random or some baseline classifier
  
  load('results/1008/POSTPILOT_HCPG_predGrp1_rounds_noo_noaddfeat.RData')
  
  # #get the standardization
  # #THIS CODE JUST FOR DOCUMENTATION; HAS BEEN DONE BEFORE AND RESULT SAVED
  # # first load postpilot data [prep for publication a new workspace]
  # setwd('C:/Users/genaucka/GitHub/PIT_GD_bv_release/R/analyses/01_classification/results/1003')
  # win_mods = cur_mod_sel_nooCV
  # win_mods = agk.recode(win_mods,c('acc'),c('ac'))
  # for (mm in 1:length(win_mods)) {
  #  pp_b_dat = featmod_coefs[[win_mods[mm]]]
  #  pp_b_dat = pp_b_dat[,grep('HCPG',names(pp_b_dat),invert=TRUE)]
  #  pp_b_dat = data.frame(pp_b_dat,pred_smoking_ftdt=dat_match$smoking_ftdt)
  #  pp_b_dat = scale(pp_b_dat)
  #  save(pp_b_dat, file=paste0('POSTPILOT_',win_mods[mm],'_stand.RData'))
  # }
  
  # predict with each model
  responses = list()
  for (mm in 1:length(list_winning_model_c_nooCV)) {
    cur_c = list_winning_model_c_nooCV[[mm]]
    cur_l = list_winning_model_l_nooCV[[mm]]
    cur_m = cur_mod_sel_nooCV[mm]
    
    mr_b_dat = featmod_coefs[[cur_m]]
    mr_b_dat = mr_b_dat[,grep('HCPG',names(mr_b_dat),invert=T)]
    
    # apply the standardization and get decision value
    load(paste0('results/1008/POSTPILOT_', cur_m,'_stand.RData'))
    pp_scale = attributes(pp_b_dat)
    
    mr_b_dat = data.frame(mr_b_dat,pred_smoking_ftdt = dat_match$smoking_ftdt)
    mr_b_dat = scale(mr_b_dat,center = pp_scale$`scaled:center`, scale = pp_scale$`scaled:scale`)
    mr_b_dat = data.frame(ones(length(mr_b_dat[,1]),1),mr_b_dat)
    
    # prediction
    responses[[mm]] = t(as.matrix(cur_c)) %*% t(as.matrix(mr_b_dat))
  }
  
  # consensus (sum of decision values, ensemble classification)
  weighted_responses = responses[[1]]
  for (mm in 2:length(responses)) {
    weighted_responses = weighted_responses + responses[[mm]]
  }
  preds     = ifelse(weighted_responses <= 0, 'HC','PG')
  
  acc = mean(preds == dat_match$HCPG)
  roc = pROC::roc(as.numeric(dat_match$HCPG),predictor=as.numeric(weighted_responses))
  auc = roc$auc
  cm  = confusionMatrix(as.factor(preds),dat_match$HCPG)
  sen = cm$byClass[1]
  spe = cm$byClass[2]
  
  # test
  p_values_ensemble_classification = list()
  p_values_ensemble_classification[['acc']] = 1-agk.density_p.c(all_accs,acc)
  p_values_ensemble_classification[['auc']] = 1-agk.density_p.c(all_aucs,auc)
  p_values_ensemble_classification[['sen']] = 1-agk.density_p.c(all_sens,sen)
  p_values_ensemble_classification[['spe']] = 1-agk.density_p.c(all_spec,spe)
  
  message('p-values for accuracy, AUC, sensitivity, specificity:')
  message(' ')
  message(p_values_ensemble_classification[['acc']])
  message(' ')
  message(p_values_ensemble_classification[['auc']])
  message(' ')
  message(p_values_ensemble_classification[['sen']])
  message(' ')
  message(p_values_ensemble_classification[['spe']])
  message(' ')
  
  # weighted models auc
  message('Applying the PIT GD classifiers to the PIT GD MRI behav data the AUC is:')
  message(' ')
  all_mod_mean_auc = auc
  message(all_mod_mean_auc)
  
  # get auc per classifier
  real_aucs = c()
  for (aa in 1:length(responses)) {
    real_aucs[aa] = as.numeric(auc(roc(as.numeric(dat_match$HCPG),as.numeric(responses[[aa]]))))
  }
  
  res = list(real_aucs=real_aucs, all_mod_mean_auc=all_mod_mean_auc, acc=acc, auc=auc, sen=sen, spe=spe,
             weighted_responses=weighted_responses, p_values_ensemble_classification=p_values_ensemble_classification)
  
  return(res)
  
}


get_aucs_of_classifiers = function() {
  
  # get the AUCs data
  
  # 1000: fMRI against control model (control is baseline model, i.e. some control variable)
  # 350: fMRI with behavioral features
  # 351: only behavioral features
  # 352: fMRI only CR
  # 353: behavioral plus fmri (only CR)
  
  aucs_env = new.env()
  aucs_env$'1000' = c('fmri_all', 'results/1000/MRI_predGrp1_rounds_wio_onlyPhys.RData', 'CV_res_list_op')
  aucs_env$'1001' = c('baseline_edu', 'results/1000/MRI_predGrp1_rounds_wio_conmod.RData', 'CVcm_res_list')
  aucs_env$'352' = c('fmri_only_cr', 'results/352/MRI_predGrp1_rounds_wio_onlyPhys.RData', 'CV_res_list_op')
  aucs_env$'350' = c('fmri_with_behav', 'results/350/MRI_predGrp1_rounds_wio_wiaddfeat.RData', 'CV_res_list')
  aucs_env$'351' = c('behav_features_only', 'results/351/MRI_predGrp1_rounds_wio_noaddfeat.RData', 'CV_res_list')
  aucs_env$'353' = c('behav_plus_fMRI_CR', 'results/353/MRI_predGrp1_rounds_wio_wiaddfeat.RData', 'CV_res_list')
  
  res = list()
  for (p in names(aucs_env)) {
    e = new.env()  # creting empty environment (i.e. dictionary)
    load(aucs_env[[p]][2], envir = e)  # loading the results
    cur_aucs = as.numeric(lapply(X=e[[aucs_env[[p]][3]]], FUN = function(x) {return(as.numeric(x$auc))}))
    res[[aucs_env[[p]][1]]] = cur_aucs
  }

  return(res)
  
}


make_df_performance_classifiers = function(list_of_aucs) {
  # the performance of the classifiers in one data-frame
  
  # param: list_of_aucs, list of aucs vectors; first one needs to be random_unbiased"
  #        following ones will be made longer if needed (NOT CUT SHORT, NOT IMPLEMENTED)
  
  # check
  stopifnot(names(list_of_aucs)[1] == 'random_unbiased')
  
  # start the data-frame
  cur_dat = data.frame(list_of_aucs$random_unbiased)
  names(cur_dat) = 'random_unbiased'
  
  for (ii in 2:length(list_of_aucs)) {
    
    # get the name of current classifier, and the data
    cur_name = names(list_of_aucs)[ii]
    cur_var_data = list_of_aucs[[ii]]
    
    # make it longer if necessary
    cur_res = agk.make.as.long(cur_dat$random_unbiased, cur_var_data)
    cur_var_data = cur_res$x2
    
    # add to data-frame
    cur_df = data.frame(cur_var_data)
    names(cur_df) = cur_name
    cur_dat = cbind(cur_dat, cur_df)
  }
  
  cur_dat              = melt(cur_dat)
  
  # SAME IDEA BUT WITH CONSTANT AS PERFORMANCE OF OLD CLASSIFIER
  #cur_dat_be = data.frame(random_classifier = all_aucs, 
  #                        mean_auc = rep(all_mod_mean_auc,length(all_aucs)),classifier = 'prev_behav_glmnet')
  
  
  #cur_dat                = rbind(cur_dat_be) #,cur_dat_gl,cur_dat_sv)
  #cur_dat                = melt(cur_dat,id.vars = c('classifier'))
  #cur_dat_H_0            = subset(cur_dat,variable == 'random_classifier')
  #cur_dat_H_0$mean_auc   = cur_dat$value[cur_dat$variable == 'mean_auc']
  #cur_dat                = cur_dat_H_0
  #cur_dat$AUC_ROC        = cur_dat$value
  #cur_dat$value          = NULL
  #cur_dat$algorithm      = cur_dat$classifier
  #cur_dat$classifier     = cur_dat$variable
  #cur_dat$classifier = agk.recode.c(cur_dat$classifier,'random_classifier','random')
  
  return(cur_dat)
}


make_density_plots = function(cur_dat) {
  
  # Plots the performance as density plots of multiple classifiers, as recorded in cur_dat
  
  
  # plot
  p = ggplot(cur_dat,aes(x=auc, y=..count../sum(..count..), fill=classifier, color=classifier)) + geom_density(alpha=0.25)
  #p = p + facet_grid(classifier ~ .) + ggtitle('AUC densities for elastic net classifier and random classifier')
  #p = p + geom_vline(aes(xintercept = as.numeric(all_mod_mean_auc)),colour = 'green',size= 1.5)
  p = p + theme_bw()
  p = p + theme(axis.text=element_text(size=14, face = "bold"),
                axis.title=element_text(size=20,face="bold"))
  p = p + coord_cartesian(xlim = c(0.4, 0.8)) 
  p = p + theme(plot.title = element_text(size=22))
  p = p + ggtitle('AUC of trained classifiers compared to random classifiers')
  p = p + xlab('AUC ROC')
  p = p + scale_y_continuous(trans="sqrt", name="probability density")
  p = p + theme(plot.title = element_text(size=25,face = 'bold'))
  p = p + theme(legend.text = element_text(size=25))
  p = p + theme(legend.title= element_text(size=25))
  print(p)
  
  
  # SAME CODE BUT WITH A CONSTANT FOR BEHAV CLASSIFIER
  #p = ggplot(cur_dat,aes(x=AUC_ROC, fill=classifier)) + geom_density(alpha=0.25)
  #p = p + facet_grid(classifier ~ .) + ggtitle('AUC densities for estimated classifier compared to random classifier')
  #p = p + ggtitle('AUC of trained classifier compared to random classifier')
  #p = p + geom_vline(aes(xintercept = mean_auc),colour = 'green',size= 1.5)
  #p = p + coord_cartesian(xlim = c(0.42, 0.8)) 
  #p = p + theme_bw()
  #p = p + theme(axis.text=element_text(size=30, face = "bold"),
  #              axis.title=element_text(size=30,face="bold"))
  #p = p + theme(plot.title = element_text(size=25,face = 'bold'))
  #p = p + theme(legend.text = element_text(size=25))
  #p = p + theme(legend.title= element_text(size=25))
  #p = p + xlab('AUC ROC')
  #print(p)
  
  return(p)
  
}




## two density plots ==============================================================






agk.make.as.long = function(x1,x2) {
  
  # utility function to repeat the shorter vector, so it reaches the lenght of the longer one
  
  if (length(x1) < length(x2)) {
    x1 = rep_len(x1,length.out = length(x2))
  } else {
    x2 = rep_len(x2,length.out = length(x1))
  }
  return(list(x1 = x1, x2 = x2))
}
