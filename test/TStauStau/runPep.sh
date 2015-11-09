#!/bin/bash

runLocalAnalysisOverSamples.py -e pepAnalysisFWLite -j tstaustau_samples_pep.json -d  /gstore/t3cms/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged/ -o /lstore/cms/cbeiraod/Pep/ -c runPep_cfg.py.templ -p "@saveSummaryTree=True @runSystematics=False @exclusiveRun=True @applyScaleFactors=True" -s 8nh

