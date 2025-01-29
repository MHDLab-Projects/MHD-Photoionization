# Convert all arguments (assumed to be SVG) to a TIFF.
# Requires Inkscape and ImageMagick 6.8 (doesn't work with 6.6.9).
# From matsen: https://gist.github.com/matsen/4263955

source ../.env

# final/figures, final/figures/SI, final/figures/schematics
TARGET_DIRS=("final/figures" "final/figures/SI" "final/figures/schematics")

for TARGET_DIR in "${TARGET_DIRS[@]}"; do
  pushd $REPO_DIR/$TARGET_DIR

  mkdir -p output

  echo "Converting SVG to PNG for $TARGET_DIR folder"

  for i in *.svg; do
    BN=$(basename $i .svg)
    inkscape --export-filename="output/$BN.png" --export-dpi 300 $i
  done

  popd
done