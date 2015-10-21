#!/bin/bash

#runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d  /gstore/t3cms/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged/ -o /lstore/cms/cbeiraod/NewObject/ -c runAnalysis_cfg.py.templ -p "@saveSummaryTree=True @runSystematics=False @exclusiveRun=True @applyScaleFactors=True @doSVfit=True" -s 8nh
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d  /gstore/t3cms/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged/ -o /lstore/cms/cbeiraod/NextToLast/ -c runAnalysis_cfg.py.templ -p "@saveSummaryTree=True @runSystematics=True @exclusiveRun=True @applyScaleFactors=True @doSVfit=False @keepOnlyPromptTaus=True" -s 8nh
runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_full.json -d  /gstore/t3cms/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged/ -o /lstore/cms/cbeiraod/NextToLast_Full/ -c runAnalysis_cfg.py.templ -p "@saveSummaryTree=True @runSystematics=True @exclusiveRun=True @applyScaleFactors=True @doSVfit=False @keepOnlyPromptTaus=False" -s 8nh

## MCClosure
#runLocalAnalysisOverSamples.py -e runTStauStauAnalysisFWLite -j tstaustau_samples_MCClosure.json -d  /gstore/t3cms/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged/ -o /lstore/cms/cbeiraod/NewObject/ -c runAnalysis_cfg.py.templ -p "@saveSummaryTree=True @runSystematics=False @exclusiveRun=True @applyScaleFactors=True @doSVfit=False" -s 8nh
