#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc472
export BUILD_ARCH=slc5_amd64_gcc462
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
cd /exper-sw/cmst3/cmssw/users/cbeiraod/SLC6/CMSSW_5_3_15/src/
eval `scramv1 runtime -sh`

cd /exper-sw/cmst3/cmssw/users/cbeiraod/SLC6/CMSSW_5_3_15/src/UserCode/llvv_fwk/test/TStauStau/Limits/
datacardMaker --json finalSelection.json --outDir ./Results/ --xsec

#Optional, comment out if unwanted
. runCombine.sh
