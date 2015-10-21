#!/bin/bash

OUTDIR=/lstore/cms/cbeiraod/DDBkgPlots
INDIR=/lstore/cms/cbeiraod/DDBkgStudy
runPlotterFWLite --noPowers --iEcm 8 --iLumi 19672 --inDir $INDIR/ --outDir $OUTDIR/ --outFile $OUTDIR/plotter.root --json tstaustau_samples_full.json --plotExt .png
# --plotExt .root

OUTDIR=/lstore/cms/cbeiraod/DDBkgFRPlots
INDIR=/lstore/cms/cbeiraod/DDBkgStudyFR
#runPlotterFWLite --noPowers --iEcm 8 --iLumi 19672 --inDir $INDIR/ --outDir $OUTDIR/ --outFile $OUTDIR/plotter.root --json tstaustau_samples_full.json --plotExt .png --plotExt .root

OUTDIR=/lstore/cms/cbeiraod/DDBkgPRPlots
INDIR=/lstore/cms/cbeiraod/DDBkgStudyPR
#runPlotterFWLite --noPowers --iEcm 8 --iLumi 19672 --inDir $INDIR/ --outDir $OUTDIR/ --outFile $OUTDIR/plotter.root --json tstaustau_samples_full.json --plotExt .png --plotExt .root

OUTDIR=/lstore/cms/cbeiraod/DDBkgPRPlotsDYOnly
INDIR=/lstore/cms/cbeiraod/DDBkgStudyPR
#runPlotterFWLite --noPowers --iEcm 8 --iLumi 19672 --inDir $INDIR/ --outDir $OUTDIR/ --outFile $OUTDIR/plotter.root --json tstaustau_samples_prompt.json --plotExt .png --plotExt .root
