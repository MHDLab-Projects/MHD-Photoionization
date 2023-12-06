: '
script to call render_scripts.sh in all subdirectories
'

render_dirs=( "absem" "tcs" "uncertainty" "various" )

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