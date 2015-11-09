#!/bin/bash

JSONFILE=tstaustau_samples_pep_plot.json

INDIR=/lstore/cms/cbeiraod/Pep
OUTDIR=/lstore/cms/cbeiraod/Pep_Plots

## 2012ABCD
runPlotterFWLite --noPowers --iEcm 8 --iLumi 19672 --inDir $INDIR/ --outDir $OUTDIR/ --outFile $OUTDIR/plotter.root --json $JSONFILE --plotExt .png --plotExt .root

