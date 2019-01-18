# ------------------------------------------
# MERGE DATA NICELY AND GET PARAMETERS FOR BEHAVIOUR
# ------------------------------------------

# patch some files together...

#NSPN
DEMO_FILE.1 = '~/Data/NSPN/demoIQ4MGM_all.csv' 
#DEMO_FILE.2 = '~/Data/NSPN/demo4BG.csv'
DEMO_FILE.2 = '~/Data/NSPN/decAc4BG2.csv'

demo.1 = read.table(DEMO_FILE.1, header=TRUE, sep=',')
demo.1$Subject = demo.1$nspnID
demo.2 = read.table(DEMO_FILE.2, header=TRUE, sep=',')
demo.2$Subject = demo.2$nspnID
demo.2$decAc = demo.2$decAc4All.bsl

#wasi_scores = c('wasi_za_vocab_raw_score', 'wasi_zl_matrix_raw_score', 'wasi_zz_iq_full2_iq')
#demo.pars.NSPN = demo.1
#inters.sub = intersect(demo.NSPN$Subject, demo.4$Subject) 
#demo.pars.NSPN[demo.NSPN$Subject %in% inters.sub,  wasi_scores] = demo.4[demo.4$Subject %in% inters.sub, wasi_scores]

demo.pars.NSPN = merge(demo.1, demo.2, by = 'Subject', all=T)

#demo.pars.NSPN$IQ = demo.pars.NSPN$wasi_zz_iq_full2_iq
#demo.pars.NSPN$IQ_vocab = demo.pars.NSPN$wasi_za_vocab_raw_score
#demo.pars.NSPN$IQ_matrix = demo.pars.NSPN$wasi_zl_matrix_raw_score

demo.pars.NSPN$age = demo.pars.NSPN$age_scan1

# transform factors to dummy
demo.pars.NSPN$sex.num = 1 * (demo.pars.NSPN$sex == 'Male')
demo.pars.NSPN$CBU = 1*(demo.pars.NSPN$centre_scan1 == 'CBU')
demo.pars.NSPN$UCL = 1*(demo.pars.NSPN$centre_scan1 == 'UCL')

# take only healthy subjects!!
#demo.pars.NSPN = subset(demo.pars.NSPN, medical_cond == 'No')

# save results
PARAMS_DEMO_FILE_NSPN = '~/Data/gng_modelling/NSPN/params_NSPN.csv'
write.table(demo.pars.NSPN, PARAMS_DEMO_FILE_NSPN, sep = ';')

# take only 2k cohort
demo.pars.2k = subset(demo.pars.NSPN, cohort = '2k_Cohort')
PARAMS_DEMO_FILE_NSPN = '~/Data/gng_modelling/NSPN/params_NSPN_2k.csv'
write.table(demo.pars.2k, PARAMS_DEMO_FILE_NSPN, sep = ';')

# those with imaging
sum(complete.cases(demo.pars.NSPN[c('Subject','age','sex','decAc')]))
