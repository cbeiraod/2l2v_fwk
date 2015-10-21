#!/bin/bash

if [[ -z "$1" ]]; then
  echo "You must specify an input directory"
  exit
fi

DIRECTORYIN1="$(readlink -m $1)"
CWD="$(pwd)"

files="$DIRECTORYIN1/*.root"

shopt -s nullglob
for file in $files; do
  echo "$file"
  root -l -b -q makeDotC.cc+\(\"$file\"\)
done
