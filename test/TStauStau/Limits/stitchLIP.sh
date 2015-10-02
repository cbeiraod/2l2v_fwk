#!/bin/bash

# StitchLIP is a script that stitches together the 3 different signal regions defined in the LIP analysis.
# The 3 signal regions are non-exclusive, but are considered in different deltaM regions

if [[ -z "$1" ]]; then
  echo "You must specify a first input directory"
  exit
fi

if [[ -z "$2" ]]; then
  echo "You must specify a second input directory"
  exit
fi

if [[ -z "$3" ]]; then
  echo "You must specify a third input directory"
  exit
fi

if [[ -z "$4" ]]; then
  echo "You must specify an output directory"
  exit
fi

DIRECTORYIN1="$(readlink -m $1)"
DIRECTORYIN2="$(readlink -m $2)"
DIRECTORYIN3="$(readlink -m $3)"
DIRECTORYOUT="$(readlink -m $4)"
CWD="$(pwd)"

if [[ -d $DIRECTORYOUT ]]; then
  echo "A results directory already exists."
  echo "Do you want to delete the directory [y/n]: (the script will stop if no is chosen)"
  read answer
  if [[ $answer == 'y' ]]; then
    echo "Deleting $DIRECTORYOUT"
    rm -Rf $DIRECTORYOUT
  else
    echo "Terminating script"
    exit
  fi
fi
mkdir $DIRECTORYOUT

for stauM in {50..500..10}
do
  for neutM in {0..480..10}
  do
    if [[ neutM -eq 0 ]]; then
      neutM=1
    fi

    if [[ neutM -ge stauM ]]; then
      continue
    fi
    let deltaM=stauM-neutM
    echo "Stau: " $stauM "; Neutralino: " $neutM "; DeltaM: " $deltaM

    file_a=""
    printf -v file_a 'SignalPoint_%d%03d.txt' "$stauM" "$neutM"
    file_b=""
    printf -v file_b 'higgsCombineS%d-N%d.Asymptotic.mH120.root' "$stauM" "$neutM"

    if [[ $deltaM -le 70 ]]; then
      cp $DIRECTORYIN1/$file_a $DIRECTORYOUT/$file_a
      cp $DIRECTORYIN1/$file_b $DIRECTORYOUT/$file_b
    else
      if [[ $deltaM -le 160 ]]; then
        cp $DIRECTORYIN2/$file_a $DIRECTORYOUT/$file_a
        cp $DIRECTORYIN2/$file_b $DIRECTORYOUT/$file_b
      else
        cp $DIRECTORYIN3/$file_a $DIRECTORYOUT/$file_a
        cp $DIRECTORYIN3/$file_b $DIRECTORYOUT/$file_b
      fi
    fi

  done
done

cd $CWD
