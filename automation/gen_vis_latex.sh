source .env

source $REPO_DIR/venv/bin/activate


if [[ "$@" =~ "data" ]]
then
    echo "Generating Dataset"
    cd $REPO_DIR/final/dataset

    for f in *.py 
    do
        python $f
    done
fi

if [[ "$@" =~ "panels" ]]
then
    echo "Generating figure panels"

    cd $REPO_DIR/final/figure_panels

    for f in *.py
    do
        python $f
    done
fi


# Figures
if [[ "$@" =~ "fig" ]]
then
    echo "Rendering figures"


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

# Collect Files
if [[ "$@" =~ "collect" ]]
then
    echo "Collecting files"

    cd $REPO_DIR/doc/paper

    cp main.pdf output/Photoionization\ Draft.pdf

    cp ../SI_man/SI_man.pdf output/Supporting\ Information.pdf

fi

cd $REPO_DIR/automation