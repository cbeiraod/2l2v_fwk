#!/bin/bash

if [[ -z "$1" ]]; then
  echo "You must specify a first input directory"
  exit
fi

if [[ -z "$2" ]]; then
  echo "You must specify a second input directory"
  exit
fi

if [[ -z "$3" ]]; then
  echo "You must specify an output directory"
  exit
fi

DIRECTORYIN1="$(readlink -m $1)"
DIRECTORYIN2="$(readlink -m $2)"
DIRECTORYOUT="$(readlink -m $3)"
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
    echo "Stau: " $stauM "; Neutralino: " $neutM

    file=""
    printf -v file 'SignalPoint_%d%03d.txt' "$stauM" "$neutM"
    file1=$DIRECTORYIN1/$file
    file2=$DIRECTORYIN2/$file

    echo $file1
    echo $file2

    option=0

    if [[ -f $file1 && -f $file2 ]]; then
      option=1
      combineCards.py $file1 $file2 &> /dev/null
      if [[ $? -ne 0 ]]; then
        combineCards.py $file1 &> /dev/null
        if [[ $? -eq 0 ]]; then
          option=2
        else
          combineCards.py $file2 &> /dev/null
          if [[ $? -eq 0 ]]; then
            option=3
          else
            option=0
          fi
        fi
      fi
    else
      if [[ -f $file1 ]]; then
        option=2
      fi
      if [[ -f $file2 ]]; then
        option=3
      fi
    fi

    case "$option" in
     1) echo " Both files exist"
        combineCards.py $file1 $file2 > $DIRECTORYOUT/$file
       ;;
     2) echo " Only the first file exists"
        combineCards.py $file1 > $DIRECTORYOUT/$file
       ;;
     3) echo " Only the second file exists"
        combineCards.py $file2 > $DIRECTORYOUT/$file
       ;;
     *) echo " No file can be processed for this point"
       ;;
    esac
  done
done

cd $CWD
