: '
script to call render_scripts.sh in all subdirectories
'

source ../.env

# Check the first argument
if [ "$1" = "final" ]; then
    expt_analysis_dir=$REPO_DIR/final
    render_dirs=( "analysis" "analysis/mwt" "analysis/topcam")
else
    expt_analysis_dir=$REPO_DIR/experiment/analysis
    render_dirs=( "mwt" "aas" "tcs" "uncertainty" "various" )
fi

purge_dirs=false
replace_existing=false

# pushd and popd are used to change directories and return to the original directory

for d in "${render_dirs[@]}"
do
    pushd $expt_analysis_dir/$d
    echo "---rendering scripts in $d---"

    ls

    if [ "$purge_dirs" = true ] ; then
        if [ -d "nb_render" ] ; then
            echo "purging nb_render directory"
            rm -rf "nb_render"
        else
            echo "nb_render does not exist"
        fi
    fi

    if [ "$replace_existing" = true ] ; then
        source $REPO_DIR/automation/render_scripts.sh -f
    else
        source $REPO_DIR/automation/render_scripts.sh
    fi

    popd
done
