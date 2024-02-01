: '
script to call render_scripts.sh in all subdirectories
'
# Default directories to process
render_dirs=( "mws" "absem" "tcs" "uncertainty" "various" )

# If an argument is provided, use it as the directories to process
if [ "$#" -gt 0 ]; then
    render_dirs=( "$@" )
fi

purge_dirs=true

for d in "${render_dirs[@]}"
do
    cd $d
    echo "---rendering scripts in $d---"

    if [ "$purge_dirs" = true ] ; then
        if [ -d "nb_render" ] ; then
            echo "purging nb_render directory"
            rm -rf "nb_render"
        else
            echo "nb_render does not exist"
        fi
    fi

    source ../render_scripts.sh
    cd ..
done