# Pipeline for main paper data and figures

auto_dir=$(pwd)
source ../.env

if [ "$MUNGE_DATA" = true ] ; then
    source munge.sh
fi

cd $auto_dir
source process_munged.sh

cd $auto_dir
source final_dataset.sh

cd $auto_dir
source gen_fig_panels.sh


cd $auto_dir
source render_figures.sh final/figures/schematics

cd $auto_dir
source render_figures.sh final/figures

