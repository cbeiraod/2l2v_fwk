#!/bin/bash

if [[ -z "$1" ]]; then
  echo "You must specify an input directory"
  exit
fi

DIRECTORY="$(readlink -m $1)"
CWD="$(pwd)"

cd $DIRECTORY

for stauM in {50..500..10}
do
  for neutM in {0..480..10}
  do
    if [[ neutM -eq 0 ]]; then
      neutM=1
    fi
    echo "Stau: " $stauM "; Neutralino: " $neutM

    file=""
    printf -v file 'SignalPoint_%d%03d.txt' "$stauM" "$neutM"
    if [[ -f $file ]]; then
      combine -M Asymptotic --run=blind $file -n S$stauM-N$neutM
    fi
  done
done

cd $CWD
