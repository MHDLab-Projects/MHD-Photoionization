# Convert all arguments (assumed to be SVG) to a TIFF.
# Requires Inkscape and ImageMagick 6.8 (doesn't work with 6.6.9).
# From matsen: https://gist.github.com/matsen/4263955

source .env

if [ -z "$1" ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

# final/figures, final/figures/SI, final/figures/schematics
TARGET_DIR=$1

pushd $REPO_DIR/$TARGET_DIR

mkdir -p output

echo "Converting SVG to PNG for $TARGET_DIR folder"

for i in *.svg; do
  BN=$(basename $i .svg)
  inkscape --export-filename="output/$BN.png" --export-dpi 300 $i
done

popd