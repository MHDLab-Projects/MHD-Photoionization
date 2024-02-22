source ../.env

source $REPO_DIR/venv/bin/activate

# Figures
if [[ "$@" =~ "fig" ]]
then
    echo "Running fig part"
    cd $REPO_DIR/final/dataset

    for f in *.py 
    do
        python $f
    done

    cd $REPO_DIR/final/figure_panels

    for f in *.py
    do
        python $f
    done

    cd $REPO_DIR/final/figures

    source gen_figs.sh
fi

# Latex 
if [[ "$@" =~ "latex" ]]
then
    echo "Running latex part"

    cd $REPO_DIR/doc/SI_man

    latexmk -synctex=1 -interaction=nonstopmode -file-line-error -pdf SI_man.tex

    source clean.sh

    cd $REPO_DIR/doc/paper

    latexmk -synctex=1 -interaction=nonstopmode -file-line-error -pdf main.tex

    source clean.sh
fi

cd $REPO_DIR/automation