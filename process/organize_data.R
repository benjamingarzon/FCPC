###############################################################################
# MERGE DATA NICELY AND GET BEHAVIOURAL PARAMETERS
###############################################################################

# patch some files together...

#NSPN
DATA_FILE.bsl = '~/Data/NSPN/demoIQ4MGM_all.csv' 
DATA_FILE.fu = '~/Data/NSPN/codes_followup.csv'
DATA_FILE.decAc = '~/Data/NSPN/decAc4b.csv' # '~/Data/NSPN/decAc4BG2.csv'
DATA_FILE.IQ = '~/Data/NSPN/wasi_IQ_Apr18.csv'


data.bsl = read.table(DATA_FILE.bsl, header=TRUE, sep=',') %>% 
  mutate(Subject = nspnID, cohort = study_primary) %>% 
  dplyr::select(Subject, 
                age_iua1, 
                cohort, 
                medical_cond, 
                medical_specify,
                centre_scan1,
                age_scan1, 
                sex) 
data.fu = read.table(DATA_FILE.fu, header=TRUE, sep='\t') %>% 
  mutate(Subject = SubjectID) %>% 
  dplyr::select(Subject,
                age_scan2, 
                centre_scan2)


# fix missing data
data.bsl[data.bsl$Subject == 21741, ]$age_scan1 = 
  data.bsl[data.bsl$Subject == 21741, ]$age_iua1

data.decAc = read.table(DATA_FILE.decAc, header=TRUE, sep=',') %>% 
  mutate(Subject = nspnID) %>% 
  rename(decAc.bsl = decAc4b.bsl, decAc.fu = decAc4b.fu2) %>% 
  dplyr::select(Subject, decAc.bsl, decAc.fu) 

# retrieve IQ values
data.IQ = read.table(DATA_FILE.IQ, header=TRUE, sep=',') %>% 
  mutate(Subject = subjectID)
data.IQ$IQcomp = as.vector(0.5*(scale(data.IQ$vocabulary_raw_score) + 
                                  scale(data.IQ$matrix_raw_score))) 

IQ.bsl = subset(data.IQ, measurement == 
                  "iua_baseline_arm_1") [c("Subject", "IQcomp")]
IQ.fu = subset(data.IQ, measurement == 
                 "iua_1st_follow_up_arm_1") [c("Subject", "IQcomp")]
IQ = merge(IQ.bsl, IQ.fu, by = "Subject", 
           all = T, suffixes = c(".bsl", ".fu"))


data.pars.NSPN = merge(data.fu, data.bsl, by = "Subject", all = T)
data.pars.NSPN = merge(data.pars.NSPN, data.decAc, by = "Subject", all = T)
data.pars.NSPN = merge(data.pars.NSPN, IQ, by = "Subject", all = T)

# fix bad naming in cohort
data.pars.NSPN$cohort[data.pars.NSPN$cohort == "2K_cohort"] = "2K_Cohort" 
data.pars.NSPN$cohort = factor(data.pars.NSPN$cohort)

# transform factors to dummy
data.pars.NSPN$sex.num = 1 * (data.pars.NSPN$sex == 'Male')
data.pars.NSPN$CBU.1 = 1*(data.pars.NSPN$centre_scan1 == 'CBU')
data.pars.NSPN$UCL.1 = 1*(data.pars.NSPN$centre_scan1 == 'UCL')
data.pars.NSPN$CBU.2 = 1*(data.pars.NSPN$centre_scan2 == 1) #different coding
data.pars.NSPN$UCL.2 = 1*(data.pars.NSPN$centre_scan2 == 2)

# take only healthy subjects!!
#data.pars.NSPN = subset(data.pars.NSPN, medical_cond == 'No')

# save results
PARAMS_DATA_FILE_NSPN = '~/Data/gng_modelling/NSPN/params_NSPN.csv'
write.table(data.pars.NSPN, PARAMS_DATA_FILE_NSPN, sep = ';')

# non-depressed - same as 2k cohort
data.pars.nondep = subset(data.pars.NSPN, cohort == '2K_Cohort') 
PARAMS_DATA_FILE_NSPN = '~/Data/NSPN/params_NSPN_nondep.csv'
write.table(data.pars.nondep, PARAMS_DATA_FILE_NSPN, sep = ';')

# create dataset with independent observations for baseline and followup
have_scan = !is.na(data.pars.nondep$age_scan1) & 
  !is.na(data.pars.nondep$age_scan2)
n_index = sum(have_scan)
index = sample(  c(rep(T, round(.4*n_index)), 
                   rep(F, n_index - round(.4*n_index))), replace = F)

data.pars.nondep.split = data.pars.nondep
data.pars.nondep.split[have_scan, c("decAc.bsl", "IQcomp.bsl", 
                                    "age_scan1")][index, ] = NA 
data.pars.nondep.split[have_scan, c("decAc.fu", "IQcomp.fu", 
                                    "age_scan2")][!index, ] = NA 

data.pars.nondep.long = rbind(
  data.pars.nondep %>% filter (!is.na(decAc.bsl)) %>% dplyr::
    select(Subject, centre_scan1, age_scan1, cohort, medical_cond, 
           medical_specify, sex, sex.num,
           decAc.bsl, IQcomp.bsl, CBU.1, UCL.1) %>% 
    rename(centre_scan = centre_scan1, age_scan = age_scan1, decAc = decAc.bsl, 
           IQcomp = IQcomp.bsl, CBU = CBU.1, UCL = UCL.1) %>% 
    mutate(group = "bsl"),
  data.pars.nondep %>% filter (!is.na(decAc.fu)) %>% 
    dplyr::select(Subject, centre_scan2, age_scan2, cohort, medical_cond, 
                  medical_specify, sex, sex.num, 
                  decAc.fu, IQcomp.fu, CBU.2, UCL.2) %>% 
    rename(centre_scan = centre_scan2, age_scan = age_scan2, decAc = decAc.fu, 
           IQcomp = IQcomp.fu, CBU = CBU.2, UCL = UCL.2) 
  %>% mutate(group = "fu")
)

PARAMS_DATA_FILE_NSPN = '~/Data/gng_modelling/NSPN/params_NSPN_nondep_long.csv'
write.table(data.pars.nondep.long, PARAMS_DATA_FILE_NSPN, sep = ';')

mydata = data.pars.nondep # data.pars.NSPN
#mydata =  data.pars.NSPN

# those with imaging
sum(complete.cases(mydata[c('Subject','age_scan1','sex','decAc.bsl')]))
sum(complete.cases(mydata[c('Subject','age_scan1','sex','decAc.bsl', 
                            'IQcomp.bsl')]))

# ica data bsl
FILE_SUBJECTS = file.path('~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica', 
                          'subjects.txt')
subjects = as.numeric(read.table(FILE_SUBJECTS))
mydata.bsl = mydata[complete.cases(mydata[c('Subject','age_scan1',
                                            'sex','decAc.bsl', 'IQcomp.bsl')]), ]
length(intersect(mydata.bsl$Subject, subjects))

# ica data fu
FILE_SUBJECTS = file.path('~/Data/NSPN/ICA_ME/nspn_ME_fu/ica200.gica', 
                          'subjects.txt')
subjects = as.numeric(read.table(FILE_SUBJECTS))
mydata.fu = mydata[complete.cases(mydata[c('Subject','age_scan2',
                                           'sex','decAc.bsl', 'IQcomp.bsl')]), ]
length(intersect(mydata.fu$Subject, subjects))
table(subset(mydata.fu, Subject %in% subjects)$cohort)
# w followup
sum(complete.cases(mydata[c('Subject','age_scan1','age_scan2',
                            'sex','decAc.bsl', 'decAc.fu')]))
sum(complete.cases(mydata[c('Subject','age_scan1','age_scan2',
                            'sex','decAc.bsl', 'decAc.fu', 'IQcomp.bsl', 
                            'IQcomp.fu')]))



# how many have followup data
#complete.fu = complete.cases(data.pars.nondep.fu[c("age_scan1", "age_scan2")])
#sum(complete.fu)
#hist(diff(t(data.pars.nondep.fu[c("age_scan1", "age_scan2")][complete.fu, ])))

#######################################
# concatenate imaging data
#######################################

create_long = function(FILE.bsl, FILE.fu, OUTFILE, sep = ","){
  
  data.bsl = read.table(FILE.bsl, sep = sep)
  data.bsl$group = "bsl"
  
  data.fu = read.table(FILE.fu, sep = sep)
  data.fu$group = "fu"
  
  data = rbind(data.bsl, data.fu)
  cc = colnames(data)[-ncol(data)]
  data = data[ c("group", cc)]
  write.table(data, file = OUTFILE, sep = sep, row.names = F, col.names = F)
  
}

if(F){
  
  create_long(
    '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica/netmats_ridge.csv',
    '~/Data/NSPN/ICA_ME/nspn_ME_fu/ica200.bsl.gica/netmats_ridge.csv',
    '~/Data/NSPN/ICA_ME/nspn_ME_long/netmats_ridge.csv'  
  )
  
  create_long(
    '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica/brain_volumes.txt',
    '~/Data/NSPN/ICA_ME/nspn_ME_fu/ica200.bsl.gica/brain_volumes.txt',
    '~/Data/NSPN/ICA_ME/nspn_ME_long/brain_volumes.txt',
    sep = " "
  )
  
  create_long(
    '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica/fd.txt',
    '~/Data/NSPN/ICA_ME/nspn_ME_fu/ica200.bsl.gica/fd.txt',
    '~/Data/NSPN/ICA_ME/nspn_ME_long/fd.txt',
    sep = " "  
  )
  
  create_long(
    '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica/dvars.txt',
    '~/Data/NSPN/ICA_ME/nspn_ME_fu/ica200.bsl.gica/dvars.txt',
    '~/Data/NSPN/ICA_ME/nspn_ME_long/dvars.txt',
    sep = " "  
  )
  
  create_long(
    '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica/meica_DOF_nom.txt',
    '~/Data/NSPN/ICA_ME/nspn_ME_fu/ica200.bsl.gica/meica_DOF_nom.txt',
    '~/Data/NSPN/ICA_ME/nspn_ME_long/meica_DOF_nom.txt',
    sep = " "  
  )
  
  create_long(
    '~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica/meica_DOF_nom.txt',
    '~/Data/NSPN/ICA_ME/nspn_ME_fu/ica200.bsl.gica/meica_DOF_nom.txt',
    '~/Data/NSPN/ICA_ME/nspn_ME_long/meica_DOF_nom.txt',
    sep = " "  
  )
  
  write.table(t(c(unlist(
    read.table('~/Data/NSPN/ICA_ME/nspn_ME/ica200.gica/subjects.txt')),
    unlist(read.table(
      '~/Data/NSPN/ICA_ME/nspn_ME_fu/ica200.bsl.gica/subjects.txt')))),
    '~/Data/NSPN/ICA_ME/nspn_ME_long/subjects.txt',
    sep = " ", row.names = F, col.names = F)  
  
}