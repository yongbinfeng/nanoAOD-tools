#!/bin/bash

INPUTDIR=/eos/cms/store/group/phys_smp/arapyan/Production_13TeV/Zmumu/ntuples/
OUTPUTDIR=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/Nano_Zmumu/
NPROC=15

RED='\033[0;31m'
NC='\033[0m'

for ifile in `ls -1 $INPUTDIR/*.root`; 
do
    if [[ "$ifile" == *"minlo"* ]]; then
        continue
    fi
    echo "*******"
    echo -e "${RED}start procesing " $ifile, $NC
    python nano_postproc.py $OUTPUTDIR $ifile --nprocesses $NPROC
    echo 
    echo
done
