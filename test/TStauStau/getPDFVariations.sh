#!/bin/bash

INDIR="/gstore/t3cms/store/user/cbeiraod/14_08_06_2l2nu_EDMtuples_merged/"
OUTDIR="/lstore/cms/cbeiraod/PDFVariations/"

cfgFile="getVariations_cfg.py.templ"

#LHAPATH="/home/cms/cbeiraod/software-area/SLC6/CMSSW_5_3_15/src/UserCode/llvv_fwk/test/TStauStau/LHAPDF/"
export LHAPATH="/exper-sw/cmst3/cmssw/users/cbeiraod/SLC6/CMSSW_5_3_15/src/UserCode/llvv_fwk/test/TStauStau/LHAPDF"

#runLocalAnalysisOverSamples.py -e computePDFvariationsFWLite -j tstaustau_samples_full.json -d ${INDIR} -o "${OUTDIR}/CT10/" -c ${cfgFile} -p "@pdfSet='CT10.LHgrid'" -s 8nh
#runLocalAnalysisOverSamples.py -e computePDFvariationsFWLite -j tstaustau_samples_full.json -d ${INDIR} -o "${OUTDIR}/MSTW/" -c ${cfgFile} -p "@pdfSet='MSTW2008nlo68cl.LHgrid'" -s 8nh
runLocalAnalysisOverSamples.py -e computePDFvariationsFWLite -j tstaustau_samples_full.json -d ${INDIR} -o "${OUTDIR}/NNPDF/" -c ${cfgFile} -p "@pdfSet='NNPDF23_nlo_as_0119.LHgrid'" -s 8nh
