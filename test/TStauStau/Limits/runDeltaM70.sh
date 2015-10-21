#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc472
export BUILD_ARCH=slc5_amd64_gcc462
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
cd /exper-sw/cmst3/cmssw/users/cbeiraod/SLC6/CMSSW_5_3_15/src/
eval `scramv1 runtime -sh`

cd /exper-sw/cmst3/cmssw/users/cbeiraod/SLC6/CMSSW_5_3_15/src/UserCode/llvv_fwk/test/TStauStau/Limits/

DIRECTORY="./TestingPlotter_Now"
#datacardMaker --verbose --json testingPlotter.json --outDir $DIRECTORY/ --xsec


DIRECTORY="./DeltaM30_Old"
#datacardMaker --verbose --json deltaM30Selection.json --outDir $DIRECTORY/ --xsec
#Optional, comment out if unwanted
#. runCombine.sh $DIRECTORY
DIRECTORY="./DeltaM70_Old"
#datacardMaker --verbose --json deltaM70Selection.json --outDir $DIRECTORY/ --xsec
#Optional, comment out if unwanted
#. runCombine.sh $DIRECTORY
DIRECTORY="./DeltaM120_Old"
#datacardMaker --verbose --json deltaM120Selection.json --outDir $DIRECTORY/ --xsec
#Optional, comment out if unwanted
#. runCombine.sh $DIRECTORY


#datacardMaker --json finalSelection.json --outDir $DIRECTORY/ --xsec

DIRECTORY="./DeltaM30_Syst"
#datacardMaker --verbose --json deltaM30Selection_Syst.json --outDir $DIRECTORY/ --xsec
#Optional, comment out if unwanted
#. runCombine.sh $DIRECTORY
DIRECTORY="./DeltaM70_Syst"
datacardMaker --verbose --json deltaM70Selection_Syst.json --outDir $DIRECTORY/ --xsec
#Optional, comment out if unwanted
. runCombine.sh $DIRECTORY
DIRECTORY="./DeltaM120_Syst"
#datacardMaker --verbose --json deltaM120Selection_Syst.json --outDir $DIRECTORY/ --xsec
#Optional, comment out if unwanted
#. runCombine.sh $DIRECTORY


DIRECTORY="./IPM_HighMT2_Old"
#datacardMaker --verbose --json IPMSelectionHighMT2.json --outDir $DIRECTORY/ --xsec
DIRECTORY="./IPM_LowMT2_Old"
#datacardMaker --verbose --json IPMSelectionLowMT2.json --outDir $DIRECTORY/ --xsec


