#!/bin/bash

OUTDIR=/lstore/cms/cbeiraod/DDBkgStudy

JSONFILE=tstaustau_samples_mc.json
INDIR=/gstore/t3cms/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged
#runLocalAnalysisOverSamples.py -e runTStauStauDDBkg -j $JSONFILE -d $INDIR/ -o $OUTDIR/ -c runAnalysis_cfg.py.templ -p "@saveSummaryTree=False @runSystematics=False @exclusiveRun=True @applyScaleFactors=False @doSVfit=False" -s 8nh


JSONFILE=tstaustau_samples_data.json
INDIR=/gstore/t3cms/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged
runLocalAnalysisOverSamples.py -e runTStauStauDDBkg -j $JSONFILE -d $INDIR/ -o $OUTDIR/ -c runAnalysis_cfg.py.templ -p "@saveSummaryTree=False @runSystematics=False @exclusiveRun=True @applyScaleFactors=False @doSVfit=False" -s 8nh
