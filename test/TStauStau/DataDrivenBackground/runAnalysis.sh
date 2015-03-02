#!/bin/bash

#Fake Rate
OUTDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/DDBkgStudyFR

JSONFILE=tstaustau_samples_full.json
INDIR=/lustre/ncg.ingrid.pt/cmst3/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged
#runLocalAnalysisOverSamples.py -e runTStauStauDDBkg -j $JSONFILE -d $INDIR/ -o $OUTDIR/ -c runAnalysis_cfg.py.templ -p "@doPrompt=False @saveSummaryTree=False @runSystematics=False @exclusiveRun=True @applyScaleFactors=False @doSVfit=False" -s 8nh


JSONFILE=tstaustau_samples_data.json
INDIR=/lustre/ncg.ingrid.pt/cmslocal/samples/tempSingleLepton
#runLocalAnalysisOverSamples.py -e runTStauStauDDBkg -j $JSONFILE -d $INDIR/ -o $OUTDIR/ -c runAnalysis_cfg.py.templ -p "@doPrompt=False @saveSummaryTree=False @runSystematics=False @exclusiveRun=True @applyScaleFactors=False @doSVfit=False" -s 8nh

#Prompt Rate
OUTDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/DDBkgStudyPR

JSONFILE=tstaustau_samples_full.json
INDIR=/lustre/ncg.ingrid.pt/cmst3/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged
runLocalAnalysisOverSamples.py -e runTStauStauDDBkg -j $JSONFILE -d $INDIR/ -o $OUTDIR/ -c runAnalysis_cfg.py.templ -p "@doPrompt=True @saveSummaryTree=False @runSystematics=False @exclusiveRun=True @applyScaleFactors=False @doSVfit=False" -s 8nh


JSONFILE=tstaustau_samples_data.json
INDIR=/lustre/ncg.ingrid.pt/cmslocal/samples/tempSingleLepton
#runLocalAnalysisOverSamples.py -e runTStauStauDDBkg -j $JSONFILE -d $INDIR/ -o $OUTDIR/ -c runAnalysis_cfg.py.templ -p "@doPrompt=True @saveSummaryTree=False @runSystematics=False @exclusiveRun=True @applyScaleFactors=False @doSVfit=False" -s 8nh
