# Custom script to clean after the paper build

echo "Cleaning paper build..."

#!/bin/sh

for ext in bbl blg log out toc bcf xml run.xml synctex.gz fls fdb_latexmk; do
    find . -type f -name "*.$ext" -delete
done