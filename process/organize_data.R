# ------------------------------------------
# MERGE DATA NICELY AND GET PARAMETERS FOR BEHAVIOUR
# ------------------------------------------

# patch some files together...
vars_to_transform = c('beta_reward', 'beta_punishment', 'alpha', 'pavlovian','noise')
vars_to_transform = c('beta_reward', 'beta_punishment', 'alpha', 'pavlovian','noise',
                      'bet_Appet', 'bet_Aver',  'lrnRate',	'PavBias',	'irNoiseXi')

myvars = c('Subject', 'bet_Appet', 'bet_Aver',  'lrnRate',	'PavBias',	'irNoiseXi',	'GoBias.x',	'iLik')

#DAD
DEMODAD_FILE = '~/Data/DAD/behaviour/demo_data.csv'
PARAMS_DAD_MCMC = '~/Data/gng_modelling/DAD/gonogo_params_DAD_mcmc.csv'
pars.DAD.MCMC = read.table(PARAMS_DAD_MCMC, sep=',')
pars.DAD.MCMC$Subject = sapply(unlist(pars.DAD.MCMC$Subject), function(x) gsub('dad', 'D', x) )
demo.DAD = read.table(DEMODAD_FILE, header=TRUE)
demo.pars.DAD = merge(demo.DAD, pars.DAD.MCMC, by = 'Subject')


#NSPN
DEMONSPN_FILE = '~/Data/NSPN/demo_data.csv'
DEMO_FILE.1 = '~/Data/NSPN/dat4MGMbsl.csv'
DEMO_FILE.2 = '~/Data/NSPN/bslMiss4MGM_edited.csv'
DEMO_FILE.3 = '~/Data/NSPN/demoIQ4MGM_all.csv' 
DEMO_FILE.4 = '~/Data/NSPN/demo4BG.csv'
DEMO_FILE.4b = '~/Data/NSPN/WASI_missing.csv'
DEMO_FILE.5 = '~/Data/NSPN/AcIn4MGM_2k.csv'
DEMO_FILE.6 = '~/Data/NSPN/decFA3BG1.csv'
DEMO_FILE.7 = '~/Data/NSPN/trancost_baseline.csv'
DEMO_FILE.8 = '~/Data/NSPN/trancost_followup.csv'
DEMO_FILE.9 = '~/Data/NSPN/codes_followup.csv'

demo.1 = read.table(DEMO_FILE.1, header=TRUE, sep=',')
demo.1$Subject = demo.1$nspnID
demo.2 = read.table(DEMO_FILE.2, header=TRUE, sep=',')
demo.2$Subject = demo.2$pt_ID
demo.3 = read.table(DEMO_FILE.3, header=TRUE, sep=',')
demo.3$Subject = demo.3$nspnID
demo.4 = read.table(DEMO_FILE.4, header=TRUE, sep=',')
demo.4$Subject = demo.4$nspnID
demo.4b = read.table(DEMO_FILE.4b, header=TRUE, sep=',')
demo.4b$Subject = demo.4b$id_nspn
demo.5 = read.table(DEMO_FILE.5, header=TRUE, sep=',')
demo.5$Subject = demo.5$ptID
demo.6 = read.table(DEMO_FILE.6, header=TRUE, sep=',')
demo.6$Subject = demo.6$nspnID
demo.7 = read.table(DEMO_FILE.7, header=TRUE, sep=',')
demo.7$Subject = demo.7$subjectID
demo.8 = read.table(DEMO_FILE.8, header=TRUE, sep=',')
demo.8$Subject = demo.8$subjectID
demo.9 = read.table(DEMO_FILE.9, header=TRUE, sep='\t')
demo.9$Subject = demo.9$SubjectID

#demo.NSPN = merge(demo.1, demo.2[, -seq(8)], by = myvars, all=T)
#demo.NSPN = merge(demo.NSPN, demo.3, by = 'Subject', all=T)

demo.NSPN = demo.3
wasi_scores = c('wasi_za_vocab_raw_score', 'wasi_zl_matrix_raw_score', 'wasi_zz_iq_full2_iq')
inters.sub = intersect(demo.NSPN$Subject, demo.4$Subject) 
demo.NSPN[demo.NSPN$Subject %in% inters.sub,  wasi_scores] = demo.4[demo.4$Subject %in% inters.sub, wasi_scores]
#demo.NSPN = merge(demo.NSPN, demo.5, by = 'Subject', all=T)
demo.NSPN = merge(demo.NSPN, demo.6, by = 'Subject', all=T)
#demo.NSPN = merge(demo.NSPN, demo.7, by = 'Subject', all=T)
#demo.NSPN = merge(demo.NSPN, demo.8, by = 'Subject', all=T)
#demo.NSPN = merge(demo.NSPN, demo.9[, c("Subject", "age_scan2")], by = 'Subject', all=T)

PARAMS_NSPN_MCMC = '~/Data/gng_modelling/NSPN/gonogo_params_NSPN_mcmc.csv'
#demo.NSPN = read.table(DEMONSPN_FILE)
pars.NSPN.MCMC = read.table(PARAMS_NSPN_MCMC, sep=',')
demo.pars.NSPN = merge(demo.NSPN, pars.NSPN.MCMC, by = 'Subject', all = T)

demo.pars.NSPN$IQ = demo.pars.NSPN$wasi_zz_iq_full2_iq
demo.pars.NSPN$IQ_vocab = demo.pars.NSPN$wasi_za_vocab_raw_score
demo.pars.NSPN$IQ_matrix = demo.pars.NSPN$wasi_zl_matrix_raw_score
demo.pars.NSPN$age = demo.pars.NSPN$age_scan1

#demo.pars.NSPN$beta_reward_t = log(demo.pars.NSPN$beta_reward)
#demo.pars.NSPN$beta_punishment_t = log(demo.pars.NSPN$beta_punishment)
#demo.pars.NSPN$alpha_t = - log(1/demo.pars.NSPN$alpha - 1)
#demo.pars.NSPN$noise_t = - log(1/demo.pars.NSPN$noise - 1)
#demo.pars.NSPN$pavlovian_t = log(demo.pars.NSPN$pavlovian)

if (F){
demo.pars.NSPN = transform_variables(demo.pars.NSPN, vars_to_transform)

# use the same transforms??
demo.pars.DAD = transform_variables(demo.pars.DAD, vars_to_transform)

# ------------------------------------------
# FIND COMPONENTS / FACTORS
# ------------------------------------------

variables.MCMC = c('beta_reward_t', 'beta_punishment_t', 'alpha_t', 'bias', 'pavlovian_t','noise_t')
variables.EM = c('bet_Appet_t', 'bet_Aver_t',  'lrnRate_t',	'GoBias.x', 'PavBias_t',	'irNoiseXi_t')

data.vars.MCMC = demo.pars.NSPN[, variables.MCMC]
#comp.MCMC = complete.cases(data.vars.MCMC)

data.vars.EM = demo.pars.NSPN[, variables.EM]
#comp.EM = complete.cases(data.vars.EM)

FA.MCMC = factanal(scale(data.vars.MCMC[comp.MCMC, ]), factors= 3, rotation='varimax', scores='regression')
FA.EM = factanal(scale(data.vars.EM[comp.EM, ]), factors= 3, rotation='varimax', scores='regression')

FA.loadings.MCMC = cor(FA.MCMC$scores, data.vars.MCMC[comp.MCMC, ])
image(FA.loadings.MCMC)
print(FA.loadings.MCMC)
FA.loadings.EM = cor(FA.EM$scores, data.vars.EM[comp.EM, ])
image(FA.loadings.EM)
print(FA.loadings.EM)

data.vars.NSPN = data.vars.MCMC
comp = comp.MCMC
FA.NSPN = FA.MCMC

#data.vars.NSPN = data.vars.EM
#comp = comp.EM
#FA.NSPN = FA.EM

pairs(FA.MCMC$loadings)
print(FA.NSPN)

library(psy)
cor(data.vars.MCMC)
cor(data.vars.EM[comp.EM, ])
#scree.plot(FA.MCMC$correlations)

  
par(mfrow=c(1,1))
#corrplot(cor(data.vars.NSPN[comp, ]), order = 'hclust', tl.col='black') 
pairs(data.vars.NSPN)


#FA.NSPN = factanal(data.vars.NSPN[comp, ], factors= 2, rotation='varimax', scores='regression')
pca.NSPN = prcomp(data.vars.NSPN[comp, ], scale.= T)
pairs(data.vars.NSPN[comp, ])
cor(data.vars.NSPN[comp, ])


demo.pars.NSPN$pca.1[comp] = pca.NSPN$x[, 1]
demo.pars.NSPN$pca.2[comp] = pca.NSPN$x[, 2]
demo.pars.NSPN$pca.3[comp] = pca.NSPN$x[, 3]
demo.pars.NSPN$pca.4[comp] = pca.NSPN$x[, 4]

demo.pars.NSPN$f.1[comp] = FA.NSPN$scores[, 1]
demo.pars.NSPN$f.2[comp] = FA.NSPN$scores[, 2]
demo.pars.NSPN$f.3[comp] = FA.NSPN$scores[, 3]

par(mfrow=c(3,1))
plot(pca.NSPN$sdev**2/sum(pca.NSPN$sdev**2), type='b')
plot(cumsum(pca.NSPN$sdev**2)/sum(pca.NSPN$sdev**2), type='b')
image(pca.NSPN$rotation)
print((pca.NSPN$rotation))

cor.test(demo.pars.NSPN$age, demo.pars.NSPN$pca.1)

} # end of finding factors


# transform factors to dummy
demo.pars.NSPN$sex.num = 1 * (demo.pars.NSPN$sex == 'Male')
demo.pars.NSPN$CBU = 1*(demo.pars.NSPN$centre_scan1 == 'CBU')
demo.pars.NSPN$UCL = 1*(demo.pars.NSPN$centre_scan1 == 'UCL')

# take only healthy subjects!!
#demo.pars.NSPN = subset(demo.pars.NSPN, medical_cond == 'No')

# save results
PARAMS_DEMO_FILE_DAD = '~/Data/gng_modelling/DAD/gonogo_params_DAD_w_demo.csv'
write.table(demo.pars.DAD, PARAMS_DEMO_FILE_DAD, sep = ';')

PARAMS_DEMO_FILE_NSPN = '~/Data/gng_modelling/NSPN/gonogo_params_NSPN_w_demo.csv'
write.table(demo.pars.NSPN, PARAMS_DEMO_FILE_NSPN, sep = ';')





