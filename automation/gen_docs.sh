source ../.env

source $REPO_DIR/venv/bin/activate

echo "Rendering main text"
cd $REPO_DIR/doc/paper
latexmk -synctex=1 -interaction=nonstopmode -file-line-error -pdf main.tex
source clean.sh

echo "Converting Main text to word document"
#TODO: where to put this?
source pandoc_convert.sh

echo "Rendering SI"
cd $REPO_DIR/doc/SI
latexmk -synctex=1 -interaction=nonstopmode -file-line-error -pdf SI.tex
source clean.sh
