#!/bin/bash

ATLAS_SUB=/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-1mm.nii.gz
ATLAS_CORT=/usr/local/fsl/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-1mm.nii.gz
MAX_SUB=`fslstats $ATLAS_SUB -R | cut -d ' ' -f2 | cut -d'.' -f2`
MAX_CORT=`fslstats $ATLAS_CORT -R | cut -d ' ' -f2 | cut -d'.' -f2`

i=200

WD=~/Data/NSPN/ICA_ME/nspn_ME/ica${i}.gica/

MELODICFILE=$WD/groupmelodic.ica/melodic_IC.nii.gz
MELODICFILE_HR=$WD/groupmelodic.ica/melodic_IC_HR.nii.gz

mri_vol2vol --mov $MELODICFILE --targ $ATLAS_SUB --regheader --o $MELODICFILE_HR
mri_segstats --id `seq $MAX_SUB` --seg $ATLAS_SUB --i $MELODICFILE_HR --avgwf $WD/groupmelodic.ica/IC_SUB.txt --excludeid 0 --empty
mri_segstats --id `seq $MAX_CORT` --seg $ATLAS_CORT --i $MELODICFILE_HR --avgwf $WD/groupmelodic.ica/IC_CORT.txt --excludeid 0 --empty

rm $MELODICFILE_HR 




