#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc472
export BUILD_ARCH=slc5_amd64_gcc462
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
cd /exper-sw/cmst3/cmssw/users/cbeiraod/SLC6/CMSSW_5_3_15/src/
eval `scramv1 runtime -sh`

cd /exper-sw/cmst3/cmssw/users/cbeiraod/SLC6/CMSSW_5_3_15/src/UserCode/llvv_fwk/test/TStauStau/Limits/


DIRECTORY="./IPM_HighMT2_Old"
#datacardMaker --verbose --json IPMSelectionHighMT2.json --outDir $DIRECTORY/ --xsec
#. runCombine.sh $DIRECTORY
DIRECTORY="./IPM_LowMT2_Old"
#datacardMaker --verbose --json IPMSelectionLowMT2.json --outDir $DIRECTORY/ --xsec
#. runCombine.sh $DIRECTORY


DIRECTORY="./IPM_HighMT2_Syst"
datacardMaker --verbose --json IPMSelectionHighMT2_Syst.json --outDir $DIRECTORY/ --xsec
. runCombine.sh $DIRECTORY
DIRECTORY="./IPM_LowMT2_Syst"
datacardMaker --verbose --json IPMSelectionLowMT2_Syst.json --outDir $DIRECTORY/ --xsec
. runCombine.sh $DIRECTORY


