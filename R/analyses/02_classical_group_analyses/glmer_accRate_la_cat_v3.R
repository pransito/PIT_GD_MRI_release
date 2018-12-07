## PREAMBLE ===================================================================
# PIT GD MRI study
# classical group-mean-differences tests/mixed ANOVAs (using lme4)
# script to describe acceptance rate
# script to fit glmer models to answer questions on acceptance rate,
# loss aversion, effects of category and group
# model comparison to get significance of overall effect of e.g. group or category
# on gain and loss
# plotting sensitivity to gain, loss and plotting loss aversion per group

# BEFORE RUNNING THIS SCRIPT
# YOU HAVE TO run the R/select_study.R script with which_study = "MRI"

## SETTINGS [nothing to change here] =============================================
# control object for the glmer fitting
cur_control = glmerControl(check.conv.grad="ignore",
                           check.conv.singular="ignore",
                           check.conv.hess="ignore",
                           optCtrl=list(optimizer = "nloptwrap",maxfun=250))
# do the boot / permutation test or just load prepared results?
# attention: it takes over an hours to run the permuation test for loss aversion group comparison
doBootPerm = 0
# bootstrap for plotting cfint's of beta gain and beta loss (already done)
doBoot     = 0
# wd for saving the results of the bootstraps
setwd(root_wd)
setwd('02_classical_group_analyses/results/effects_under_0_la')
if (which_study == 'MRI') {
  setwd('MRI')
} else if (which_study == 'POSTPILOT_HCPG') {
  setwd('behav')
} else {
  stop('which_study is set to an unknown value!')
}
bootResWd = getwd()
# how many bootstraps/permutations?
cur_num   = 300
# how many cpus to use?
cur_cpus  = detectCores()-1
# fit glmer models?
# careful this takes on an intel core i7 with 8GB RAM about 10 minutes
doFitGlmer = 1
# put in the original fixed effect estimation
put_in_original_fe = 1

## glmer models la ============================================================
if (doFitGlmer) {
  modla_00  = glmer(accept_reject ~ 1 + (1|subject) + (1|stim) + (1|cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_01  = glmer(accept_reject ~ gain + loss  + (gain + loss |subject) + (gain + loss |stim) + (gain + loss |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_0g  = glmer(accept_reject ~ (gain + loss )*HCPG + (gain + loss |subject) + (gain + loss |stim)  + (gain + loss |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_c0  = glmer(accept_reject ~ (gain + loss ) + cat + (gain + loss + cat|subject) + (gain + loss |stim) + (gain + loss |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_cg  = glmer(accept_reject ~ (gain + loss )*HCPG + cat*HCPG + (gain + loss + cat |subject) + (gain + loss |stim)  + (gain + loss |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_cgi = glmer(accept_reject ~ (gain + loss )*cat*HCPG + ((gain + loss )*cat|subject) + (gain + loss |stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_ci  = glmer(accept_reject ~ (gain + loss )*cat + ((gain + loss )*cat|subject) + (gain + loss |stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
}

## glmer models lae (MRI study) ===============================================
if (doFitGlmer) {
  modlae_00  = glmer(accept_reject ~ (gain + loss + ed_abs) + (gain + loss + ed_abs|subject) + (gain + loss + ed_abs |stim) + (gain + loss + ed_abs |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modlae_0g  = glmer(accept_reject ~ (gain + loss + ed_abs)*HCPG + (gain + loss + ed_abs|subject) + (gain + loss + ed_abs |stim) + (gain + loss + ed_abs |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  
  modlae_c0  = glmer(accept_reject ~ (gain + loss + ed_abs) + cat + (gain + loss + ed_abs + cat|subject) + (gain + loss + ed_abs |stim) + (gain + loss + ed_abs |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modlae_cg  = glmer(accept_reject ~ (gain + loss + ed_abs)*HCPG + cat*HCPG + (gain + loss + ed_abs + cat |subject) + (gain + loss + ed_abs |stim)  + (gain + loss + ed_abs |cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
}

## check model fit per subject
cur_dp         = modla_cg@frame
cur_dp$pred_00 = as.numeric(as.numeric(predict(modla_00) >= 0) == cur_dp$accept_reject)
cur_dp$pred_cg = as.numeric(as.numeric(predict(modla_cg) >= 0) == cur_dp$accept_reject)
dens_df        = aggregate(cbind(pred_00,pred_cg) ~ subject + HCPG, data = cur_dp, FUN = mean)
message('The model modla_cg makes better predictions than modla_00 in ...')
message(' ')
message(mean(dens_df$pred_cg > dens_df$pred_00)*100)
message(' ')
message('... percent subjects.')

# acceptance rate based on laec_group model (MRI study) #######################
# general acceptance rate difference between group
laecg_effects = summary(modlae_cg)
#laecg_effects = summary(modla_0g)

HC_acc   = logistic(laecg_effects$coefficients['(Intercept)',1])
HC_acc_u = logistic(laecg_effects$coefficients['(Intercept)',1] + 1.96*laecg_effects$coefficients['(Intercept)',2])
HC_acc_l = logistic(laecg_effects$coefficients['(Intercept)',1] - 1.96*laecg_effects$coefficients['(Intercept)',2])

GD_acc   = logistic(laecg_effects$coefficients['(Intercept)',1] + laecg_effects$coefficients['HCPGPG',1])
GD_acc_u = logistic(laecg_effects$coefficients['(Intercept)',1] + laecg_effects$coefficients['HCPGPG',1] + 1.96*laecg_effects$coefficients['HCPGPG',2])
GD_acc_l = logistic(laecg_effects$coefficients['(Intercept)',1] + laecg_effects$coefficients['HCPGPG',1] - 1.96*laecg_effects$coefficients['HCPGPG',2])

message('Acceptance rate for HC, for GD at neutral category and given all other predictors are 0:')
message(HC_acc)
message(GD_acc)
HC_gain = laecg_effects$coefficients['gain',1]
GD_gain = HC_gain + laecg_effects$coefficients['gain:HCPGPG',1]
HC_loss = laecg_effects$coefficients['loss',1]
GD_loss = HC_loss + laecg_effects$coefficients['loss:HCPGPG',1]

# sensitivity to gain, loss, LA
message('Sensitivity to gain, loss for HC, GD was: ')
message('gain')
message(HC_gain)
message(GD_gain)
message('loss')
message(HC_loss)
message(GD_loss)
message('Loss aversion for HC, GD was according to fixed effects of gain and loss: ')
message(-HC_loss/HC_gain)
message(-GD_loss/GD_gain)

# acceptance rate per category controlled for neutral
HC_acc_gam = logistic(laecg_effects$coefficients['catgambling',1] + laecg_effects$coefficients['(Intercept)',1]) - HC_acc
HC_acc_neg = logistic(laecg_effects$coefficients['catnegative',1] + laecg_effects$coefficients['(Intercept)',1]) - HC_acc 
HC_acc_pos = logistic(laecg_effects$coefficients['catpositive',1] + laecg_effects$coefficients['(Intercept)',1]) - HC_acc

HC_acc_gam_u = logistic(laecg_effects$coefficients['catgambling',1] + laecg_effects$coefficients['(Intercept)',1] + 1.96*laecg_effects$coefficients['catgambling',2]) - HC_acc
HC_acc_neg_u = logistic(laecg_effects$coefficients['catnegative',1] + laecg_effects$coefficients['(Intercept)',1] + 1.96*laecg_effects$coefficients['catnegative',2]) - HC_acc 
HC_acc_pos_u = logistic(laecg_effects$coefficients['catpositive',1] + laecg_effects$coefficients['(Intercept)',1] + 1.96*laecg_effects$coefficients['catpositive',2]) - HC_acc

HC_acc_gam_l = logistic(laecg_effects$coefficients['catgambling',1] + laecg_effects$coefficients['(Intercept)',1] - 1.96*laecg_effects$coefficients['catgambling',2]) - HC_acc
HC_acc_neg_l = logistic(laecg_effects$coefficients['catnegative',1] + laecg_effects$coefficients['(Intercept)',1] - 1.96*laecg_effects$coefficients['catnegative',2]) - HC_acc 
HC_acc_pos_l = logistic(laecg_effects$coefficients['catpositive',1] + laecg_effects$coefficients['(Intercept)',1] - 1.96*laecg_effects$coefficients['catpositive',2]) - HC_acc

GD_acc_gam = logistic(laecg_effects$coefficients['catgambling',1] + laecg_effects$coefficients['(Intercept)',1] +laecg_effects$coefficients['HCPGPG:catgambling',1]) - GD_acc
GD_acc_neg = logistic(laecg_effects$coefficients['catnegative',1] + laecg_effects$coefficients['(Intercept)',1] +laecg_effects$coefficients['HCPGPG:catnegative',1]) - GD_acc
GD_acc_pos = logistic(laecg_effects$coefficients['catpositive',1] + laecg_effects$coefficients['(Intercept)',1] +laecg_effects$coefficients['HCPGPG:catpositive',1]) - GD_acc

GD_acc_gam_u = logistic(laecg_effects$coefficients['catgambling',1] + laecg_effects$coefficients['(Intercept)',1] +laecg_effects$coefficients['HCPGPG:catgambling',1] + 1.96*laecg_effects$coefficients['HCPGPG:catgambling',2]) - GD_acc
GD_acc_neg_u = logistic(laecg_effects$coefficients['catnegative',1] + laecg_effects$coefficients['(Intercept)',1] +laecg_effects$coefficients['HCPGPG:catnegative',1] + 1.96*laecg_effects$coefficients['HCPGPG:catnegative',2]) - GD_acc
GD_acc_pos_u = logistic(laecg_effects$coefficients['catpositive',1] + laecg_effects$coefficients['(Intercept)',1] +laecg_effects$coefficients['HCPGPG:catpositive',1] + 1.96*laecg_effects$coefficients['HCPGPG:catpositive',2]) - GD_acc

GD_acc_gam_l = logistic(laecg_effects$coefficients['catgambling',1] + laecg_effects$coefficients['(Intercept)',1] +laecg_effects$coefficients['HCPGPG:catgambling',1] - 1.96*laecg_effects$coefficients['HCPGPG:catgambling',2]) - GD_acc
GD_acc_neg_l = logistic(laecg_effects$coefficients['catnegative',1] + laecg_effects$coefficients['(Intercept)',1] +laecg_effects$coefficients['HCPGPG:catnegative',1] - 1.96*laecg_effects$coefficients['HCPGPG:catnegative',2]) - GD_acc
GD_acc_pos_l = logistic(laecg_effects$coefficients['catpositive',1] + laecg_effects$coefficients['(Intercept)',1] +laecg_effects$coefficients['HCPGPG:catpositive',1] - 1.96*laecg_effects$coefficients['HCPGPG:catpositive',2]) - GD_acc

# make a graph of this
HC_gam = c(HC_acc_gam,HC_acc_gam_l,HC_acc_gam_u,'gambling','HC')
HC_neg = c(HC_acc_neg,HC_acc_neg_l,HC_acc_neg_u,'negative','HC')
HC_pos = c(HC_acc_pos,HC_acc_pos_l,HC_acc_pos_u,'positive','HC')

GD_gam = c(GD_acc_gam,GD_acc_gam_l,GD_acc_gam_u,'gambling','GD')
GD_neg = c(GD_acc_neg,GD_acc_neg_l,GD_acc_neg_u,'negative','GD')
GD_pos = c(GD_acc_pos,GD_acc_pos_l,GD_acc_pos_u,'positive','GD')

cur_mat = as.data.frame(rbind(HC_gam,HC_neg,HC_pos,GD_gam,GD_neg,GD_pos))
names(cur_mat) = c('delta_percentage','CI_lower','CI_upper','category','group')
cur_fun= function(x) {return(as.numeric(as.character(x)))}
cur_mat[1:3]  = lapply(cur_mat[1:3],FUN = cur_fun)

mRat  = ggplot(cur_mat, aes(category, delta_percentage,fill=group))
mRat  = mRat + labs(x='category', y=paste('shift in acceptance rate (',0.95*100,'% CI)'))
mRat  = mRat + ggtitle("Shift in mean acceptance across categories")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), position=dodge, width=0.25) + theme_bw()
mRat  = mRat + theme(axis.text=element_text(size=18),
                     axis.title=element_text(size=20,face="bold")) + theme(plot.title = element_text(size=22)) 
mRat  = mRat +  theme(legend.title = element_text(size=18)) +  theme(legend.text = element_text(size=15))
print(mRat)

## loss aversion (la) overall and group comparison ============================
# stats glmer
print(anova(modla_00,modla_01,modla_0g,modla_cg,modla_cgi)) # all models
print(anova(modla_00,modla_01,modla_0g))                    # only loss aversion relevant
print(anova(modlae_c0,modlae_cg))                           # simple comparison group using the laec model; for MRI paper
print(anova(modlae_00,modlae_0g))                           # simple comparison group using the lae model; for MRI paper
print(summary(modlae_cg))