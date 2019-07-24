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
# YOU HAVE TO run the R/select_study.R script with which_study = "POSTPILOT_HCPG"

## SETTINGS [nothing to change here ==============================================
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
cur_num   = 100
# how many cpus to use?
cur_cpus  = detectCores()-1
# fit glmer models?
# careful this takes on an intel core i7 with 8GB RAM about 10 minutes
doFitGlmer = 1
# put in the original fixed effect estimation
put_in_original_fe = 1

# save the df
data_pdt_cat = data_pdt

## glmer models acceptance rate ===============================================
if (doFitGlmer) {
  moda_00  = glmer(accept_reject ~ 1 + (1|subject) + (1|stim) + (1|cat),data = data_pdt,family = 'binomial')
  moda_01  = glmer(accept_reject ~ cat + (cat|subject) + (1|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  moda_02  = glmer(accept_reject ~ cat*HCPG + (cat|subject) + (1|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  moda_01b = glmer(accept_reject ~ HCPG + (1|subject) + (1|stim) + (1|cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
}

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
if (doFitGlmer & which_study == 'MRI') {
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

## acceptance rate and under different cue conditions #########################
# acceptance rate graph, descriptives (CIs over subjects per category)
# acceptance rate not based on lmer model
mod_acc        = aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject,data_pdt$cat), FUN=mean.rmna)
names(mod_acc) = c('subject','category','mean_acceptance')
mod_acc$Group  = agk.recode.c(mod_acc$subject,dat_match$VPPG,dat_match$HCPG)
mod_acc        = aggregate(mod_acc$mean_acceptance,by=list(mod_acc$Group,mod_acc$cat),FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
mod_acc        = data.frame(mod_acc[[1]], mod_acc[[2]],mod_acc[[3]])
names(mod_acc) = c('Group','category','mean_acceptance','ci_0025','ci_0975')
mod_acc$Group  = agk.recode.c(mod_acc$Group,'PG','GD')

mRat  = ggplot(mod_acc, aes(category, mean_acceptance,fill=Group))
mRat  = mRat + labs(x='category', y=paste('Mean of acceptance (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Mean acceptance across categories")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = ci_0025, ymax = ci_0975), position=dodge, width=0.25) + theme_bw()
mRat  = agk.make.plot.text.bigger(mRat)
print(mRat)

# effect of cat*HCPG
message('Effect of category X GROUP')
print(anova(moda_01,moda_02))
print(summary(moda_02))

## acceptance rate and under different cue conditions indivd. #################
mod_acc             = aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject,data_pdt$cat), FUN=mean.rmna)
names(mod_acc)      = c('subject','category','mean_acceptance')
mod_acc$Group       = agk.recode.c(mod_acc$subject,dat_match$VPPG,dat_match$HCPG)
mod_acc$Group       = agk.recode.c(mod_acc$Group,'PG','GD')
mod_acc_agg         = aggregate(mod_acc$mean_acceptance,by=list(mod_acc$Group,mod_acc$cat),FUN=mean)
names(mod_acc_agg)  = c('Group','category','mean_acceptance')

mRat  = ggplot(mod_acc, aes(category, mean_acceptance,color=Group)) 
mRat  = mRat + geom_line(aes(group = subject), alpha = .3)
mRat  = mRat + labs(x='category', y=paste('Mean of acceptance'))
mRat  = mRat + geom_line(data = mod_acc_agg, alpha = .9, size = 3, aes(group = Group))
mRat  = mRat + ggtitle("Mean acceptance across categories:\nindividual data and grand means by group")
mRat  = mRat +  theme_bw()
mRat  = agk.make.plot.text.bigger(mRat)
mRat

# acceptance rate based on laec_group model (MRI study) #######################
# general acceptance rate difference between group (using the neutral cat)
if (which_study == 'MRI') {
  message('acceptance rate based on laec_group model (MRI study)')
  laecg_effects = summary(modlae_cg)

  HC_acc   = logistic(laecg_effects$coefficients['(Intercept)',1])
  HC_acc_u = logistic(laecg_effects$coefficients['(Intercept)',1] + 1.96*laecg_effects$coefficients['(Intercept)',2])
  HC_acc_l = logistic(laecg_effects$coefficients['(Intercept)',1] - 1.96*laecg_effects$coefficients['(Intercept)',2])
  
  GD_acc   = logistic(laecg_effects$coefficients['(Intercept)',1] + laecg_effects$coefficients['HCPGPG',1])
  GD_acc_u = logistic(laecg_effects$coefficients['(Intercept)',1] + laecg_effects$coefficients['HCPGPG',1] + 1.96*laecg_effects$coefficients['HCPGPG',2])
  GD_acc_l = logistic(laecg_effects$coefficients['(Intercept)',1] + laecg_effects$coefficients['HCPGPG',1] - 1.96*laecg_effects$coefficients['HCPGPG',2])
  
  message('Acceptance rate for HC, for GD at neutral category and given all other predictors are 0:')
  message(HC_acc)
  message(GD_acc)
  
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
  
  # print delta ac rate for GD and HC
  message('Acceptance rate delta for HC, for GD at gambling category (vs. neutral):')
  message(HC_acc_gam)
  message(GD_acc_gam)
  message('Acceptance rate delta for HC, for GD at negative category (vs. neutral):')
  message(HC_acc_neg)
  message(GD_acc_neg)
  message('Acceptance rate delta for HC, for GD at positive category (vs. neutral):')
  message(HC_acc_pos)
  message(GD_acc_pos)
  
  
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
  #cur_mat = melt(cur_mat,id.vars = c('category','group'))
  
  mRat  = ggplot(cur_mat, aes(category, delta_percentage,fill=group))
  mRat  = mRat + labs(x='category', y=paste('shift in acceptance rate (',0.95*100,'% CI)'))
  mRat  = mRat + ggtitle("Shift in mean acceptance across categories\n(relative to neutral, based on laecg model)")
  mRat  = mRat + geom_bar(position="dodge", stat="identity")
  dodge = position_dodge(width=0.9)
  mRat  = mRat + geom_bar(position=dodge, stat="identity")
  mRat  = mRat + geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), position=dodge, width=0.25) + theme_bw()
  mRat  = agk.make.plot.text.bigger(mRat)
  print(mRat)
  
  print(anova(modlae_c0,modlae_cg))
  print(summary(modlae_cg))
}

## sensitivity to gain, loss, LA ==============================================
laecg_effects = summary(modla_0g)
  
HC_gain = laecg_effects$coefficients['gain',1]
GD_gain = HC_gain + laecg_effects$coefficients['gain:HCPGPG',1]
HC_loss = laecg_effects$coefficients['loss',1]
GD_loss = HC_loss + laecg_effects$coefficients['loss:HCPGPG',1]  

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

print (anova(modla_00,modla_01,modla_0g,modla_cg,modla_cgi))
print(summary(modla_cg))

if (which_study == 'MRI') {
  # LA effects based on model with ED
  print (anova(modlae_00,modlae_0g))
}

if (which_study == 'MRI') {
  # LA effects overall; no group effect
  laecg_effects = summary(modla_01)
  
  HC_gain = laecg_effects$coefficients['gain',1]
  HC_loss = laecg_effects$coefficients['loss',1]
  
  message('Sensitivity to gain, loss for all subs: ')
  message('gain')
  message(HC_gain)
  message('loss')
  message(HC_loss)
  message('Loss aversion for all subs according to fixed effects of gain and loss: ')
  message(-HC_loss/HC_gain)
}

if (which_study == 'MRI') {
  # LA effects overall; no group effect
  laecg_effects = summary(modlae_cg)
  
  HC_gain = laecg_effects$coefficients['gain',1]
  GD_gain = HC_gain + laecg_effects$coefficients['gain:HCPGPG',1]
  HC_loss = laecg_effects$coefficients['loss',1]
  GD_loss = HC_loss + laecg_effects$coefficients['loss:HCPGPG',1]  
  
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
}

# acceptance rate only between group
mod_accnc        = aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject), FUN=mean.rmna)
names(mod_accnc) = c('subject','mean_acceptance')
mod_accnc$Group  = agk.recode.c(mod_accnc$subject,dat_match$VPPG,dat_match$HCPG)
mod_accnc        = aggregate(mod_accnc$mean_acceptance,by=list(mod_accnc$Group),FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
mod_accnc        = data.frame(mod_accnc[[1]],mod_accnc[[2]])
names(mod_accnc) = c('Group','mean_acceptance','ci_0025','ci_0975')
mod_accnc$Group  = agk.recode.c(mod_accnc$Group,'PG','GD')
message('Acceptance rate per group descriptively across all categories.')
print(mod_accnc)

message('effect of category and group')
print(anova(moda_00,moda_01,moda_02)) # effect of category and group

message('stats without cat (simple acceptance rate difference between groups)')
print(anova(moda_00,moda_01b))

## get the distribution of KFG and SOGS
KFG_SOGS = dat_match[c('HCPG','KFG','SOGS')]
KFG_SOGS[c('KFG','SOGS')] = KFG_SOGS[c('KFG','SOGS')] + randn(dim(KFG_SOGS)[1],dim(KFG_SOGS)[2]-1)*0.1
KFG_SOGS$Group= as.factor(agk.recode(as.character(KFG_SOGS$HCPG),c('HC','PG'),c('HC','GD')))

p = ggplot(KFG_SOGS, aes(KFG, SOGS))
p = p + geom_point(size = 4, shape = 17)
p = p + geom_point(aes(colour = Group), size = 4, shape = 17)
p = p + xlab('Score Short Questionnaire Pathological Gambling (KFG)')
p = p + ylab('Score South Oaks Gambling Screen (SOGS)')
p = p + geom_hline(yintercept = 4) # linetype, color, size)
p = p + geom_vline(xintercept = 15) # linetype, color, size)
p = p + ggtitle('Sorting subjects into groups according to score on the KFG')
p = p + theme(axis.text=element_text(size=14, face = "bold"),
              axis.title=element_text(size=20,face="bold"))
p = p + theme(plot.title = element_text(size=22))
p = p + theme(legend.text = element_text(size=18))
p = p + theme(legend.title= element_text(size=18))
print(p)
## get cfint using bootmer ====================================================
# does not work or is wayyyy too slow
#confint.merMod(modla_0g,method = 'Wald')
