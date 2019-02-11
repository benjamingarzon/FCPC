# ------------------------------------------
# MERGE DATA NICELY AND GET PARAMETERS FOR BEHAVIOUR
# ------------------------------------------

# patch some files together...

#NSPN
DEMO_FILE.1 = '~/Data/NSPN/demoIQ4MGM_all.csv' 
DEMO_FILE.2 = '~/Data/NSPN/decAc4BG2.csv' # '~/Data/NSPN/demo4BG.csv'

demo.1 = read.table(DEMO_FILE.1, header=TRUE, sep=',')
demo.1$Subject = demo.1$nspnID
demo.2 = read.table(DEMO_FILE.2, header=TRUE, sep=',')
demo.2$Subject = demo.2$nspnID
demo.2$decAc = demo.2$decAc4All.bsl

# retrieve IQ values
wasi_scores = c('wasi_za_vocab_raw_score', 'wasi_zl_matrix_raw_score', 'wasi_zz_iq_full2_iq')
DEMO_FILE.3 = '~/Data/NSPN/demoIQ4MGM_all.csv' 
DEMO_FILE.4 = '~/Data/NSPN/demo4BG.csv'
DEMO_FILE.5 = '~/Data/NSPN/WASI_missing.csv'
demo.3 = read.table(DEMO_FILE.3, header=TRUE, sep=',') %>% mutate(Subject = nspnID) 
demo.3 = demo.3[c("Subject", wasi_scores)]
demo.4 = read.table(DEMO_FILE.4, header=TRUE, sep=',') %>% mutate(Subject = nspnID)
demo.4 = demo.4[c("Subject", wasi_scores)]
demo.5 = read.table(DEMO_FILE.5, header=TRUE, sep=',') %>% mutate(Subject = id_nspn)
demo.5 = demo.5[c("Subject", wasi_scores)]

# join and reconcile conflicts
demo.3 = subset(demo.3, !Subject %in% c(19935, 27813, 35501))# these seem incorrect in the first file
demo.IQ = rbind(demo.3, demo.4, demo.5)
demo.IQ = distinct(demo.IQ[ !is.na(rowSums(demo.IQ[wasi_scores])), ])

#subset(demo.IQ, Subject %in% names(which(table(demo.IQ$Subject)>1)))

demo.pars.NSPN = merge(demo.1, demo.2, by = 'Subject', all=T)
demo.pars.NSPN = merge(demo.pars.NSPN %>% dplyr::select(-one_of(wasi_scores)) , demo.IQ, by = 'Subject', all=T)

demo.pars.NSPN = demo.pars.NSPN %>% mutate(IQ = wasi_zz_iq_full2_iq,
                                           IQvocab = wasi_za_vocab_raw_score,
                                           IQmatrix = wasi_zl_matrix_raw_score,
                                           age = age_scan1)

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
sum(complete.cases(demo.pars.NSPN[c('Subject','age','sex','decAc', wasi_scores)]))
