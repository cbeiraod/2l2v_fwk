#!/bin/bash

if [[ -z "$1" ]]; then
  echo "You must specify a json file"
  exit
fi

FILE="$(readlink -m $1)"

TEMP="grep 'iLumi' $FILE | wc -l"
NRounds="$(eval $TEMP)"
CWD="$(pwd)"
CWD="$(readlink -m $CWD)"

OUTDIR="$HOME/local-area/NewCutOptim/Step5"
OUTDIR="$HOME/local-area/NewCutOptim/FullSimplifiedFinal"
OUTDIR="$HOME/local-area/NewCutOptim/NMinus1_take2"
OUTDIR="$HOME/local-area/NewCutOptim/FullSimplifiedFinal_Syst"
#OUTDIR="$HOME/local-area/NewCutOptim/Full1"
OUTDIR="$(readlink -m $OUTDIR)"

echo "Found $NRounds rounds in file $1"
if [[ $NRounds == 0 ]]; then
  echo "There are no rounds to submit jobs for"
  exit
fi

if [[ -d $OUTDIR ]]; then
  echo "A results directory already exists."
  echo "Do you want to delete the directory [y/n]: (the script will stop if no is chosen)"
  read answer
  if [[ $answer == 'y' ]]; then
    echo "Deleting $OUTDIR"
    rm -Rf $OUTDIR
  else
    echo "Terminating script"
    exit
  fi
fi
mkdir $OUTDIR
TEMP="cp $FILE $OUTDIR/cutOptim.json"
eval $TEMP

for (( round=0; round<$NRounds; round++ ))
do
  echo "Processing Round " $round

  TEMP="cp runRound.sh.templ $OUTDIR/runRound$round.sh"
  eval $TEMP

  sed -i -e "s|#ROUND|$round|g" "$OUTDIR/runRound$round.sh"
  sed -i -e "s|#CWD|$CWD|g" "$OUTDIR/runRound$round.sh"
  sed -i -e "s|#FILE|$FILE|g" "$OUTDIR/runRound$round.sh"
  sed -i -e "s|#OUTDIR|$OUTDIR|g" "$OUTDIR/runRound$round.sh"
done

echo "Submitting the jobs"
for (( round=0; round<$NRounds; round++ ))
do
  TEMP="qsub $OUTDIR/runRound$round.sh"
  eval $TEMP
done
