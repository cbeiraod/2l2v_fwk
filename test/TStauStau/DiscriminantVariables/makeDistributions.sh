#!/bin/bash

OUTDIR=/lstore/cms/cbeiraod/Distributions
INDIR=/lstore/cms/cbeiraod/NewObject
OUTDIR=./

INDIR=/lstore/cms/cbeiraod/FinalUpdated_Full
#makeDistributions --json tstaustau_samples_MC.json --variables variables.json --inDir $INDIR/ --outDir $OUTDIR/ --signalSelection "d_stauMass-d_neutralinoMass==50" --extraSignal "d_stauMass-d_neutralinoMass==200" --plotExt .png --plotExt .C --baseSelection "b_selected" --iLumi 19672 --pointVar "d_stauMass*1000+d_neutralinoMass" --sigXSec 0.1 --sigNInitEvents 10000 --unblind --weightVar "d_weight" --puweightVar "d_PUweight"

INDIR=/lstore/cms/cbeiraod/FinalUpdated
makeDistributions --json tstaustau_samples_full.json --variables variables.json --inDir $INDIR/ --outDir $OUTDIR/ --signalSelection "d_stauMass-d_neutralinoMass==50" --extraSignal "d_stauMass-d_neutralinoMass==200" --plotExt .png --plotExt .C --baseSelection "b_selected" --iLumi 19672 --pointVar "d_stauMass*1000+d_neutralinoMass" --sigXSec 0.1 --sigNInitEvents 10000 --unblind --weightVar "d_weight" --puweightVar "d_PUweight"
