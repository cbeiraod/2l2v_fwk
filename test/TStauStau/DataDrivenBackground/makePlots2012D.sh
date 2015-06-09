#!/bin/bash

OUTDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/DDBkgFRPlots
INDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/DDBkgStudyFR
runPlotterFWLite --noPowers --iEcm 8 --iLumi 7329 --inDir $INDIR/ --outDir $OUTDIR/ --outFile $OUTDIR/plotter.root --json tstaustau_samples_2012D.json --plotExt .png --plotExt .root

OUTDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/DDBkgPRPlots
INDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/DDBkgStudyPR
#runPlotterFWLite --noPowers --iEcm 8 --iLumi 7329 --inDir $INDIR/ --outDir $OUTDIR/ --outFile $OUTDIR/plotter.root --json tstaustau_samples_2012D.json --plotExt .png --plotExt .root

OUTDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/DDBkgPRPlotsDYOnly
INDIR=/lustre/ncg.ingrid.pt/cmslocal/cbeiraod/DDBkgStudyPR
#runPlotterFWLite --noPowers --iEcm 8 --iLumi 19672 --inDir $INDIR/ --outDir $OUTDIR/ --outFile $OUTDIR/plotter.root --json tstaustau_samples_prompt.json --plotExt .png --plotExt .root
