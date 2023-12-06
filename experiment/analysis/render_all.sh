: '
script to call render_scripts.sh in all subdirectories
'

render_dirs=( "absem" "tcs" "uncertainty" "various" )

for d in "${render_dirs[@]}"
do
    cd $d
    echo "---rendering scripts in $d---"
    source ../render_scripts.sh
    cd ..
done