# deltes data directories, except for munged data

# experiment\data\proc_data
# final\dataset\output
# final\dataset\topcam\output
# modeling\viability\dataset\output

# final\figure_panels\output\figures

# iterate through all folders and remove contents

dir_list=('experiment/data/proc_data' 'final/dataset/output' 'final/dataset/topcam/output' 'modeling/viability/dataset/output' 'final/figure_panels/output/figures')

for dir in ${dir_list[@]}
do
    rm -r $REPO_DIR/$dir/*
done

