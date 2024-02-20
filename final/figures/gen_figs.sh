#!/bin/sh
# Convert all arguments (assumed to be SVG) to a TIFF.
# Requires Inkscape and ImageMagick 6.8 (doesn't work with 6.6.9).
# From matsen: https://gist.github.com/matsen/4263955

mkdir -p output

echo "Converting SVG to PNG for figures folder"

for i in *.svg; do
  BN=$(basename $i .svg)
  inkscape --export-filename="output/$BN.png" --export-dpi 300 $i
done

cd SI

mkdir -p output

echo "Converting SVG to PNG for SI folder"


for i in *.svg; do
  echo $i
  BN=$(basename $i .svg)
  inkscape --export-filename="output/$BN.png" --export-dpi 300 $i
done

cd ..