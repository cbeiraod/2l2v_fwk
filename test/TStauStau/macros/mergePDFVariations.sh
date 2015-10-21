#!/bin/bash

cwd=$(pwd)

cd ~/local-area/PDFVariations/
for file in CT10/*.root;
do
  currentFile=$(basename $file)
  echo "Merging $currentFile"
  hadd $currentFile CT10/$currentFile MSTW/$currentFile NNPDF/$currentFile
done

cd $cwd
