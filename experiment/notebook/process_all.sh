#
!/bin/bash
``
source ../../.env

nb_dir=$EXPT_NB_DIR

dates=("2023-04-07" "2023-05-12" "2023-05-18" "2023-05-24")

process_mws=true
process_absem=true

purge_dirs=false

for date in ${dates[*]};

do
    echo --------$date---------

    date_dir=$nb_dir/$date

    if $process_mws; then

        cd $date_dir/mws

        if $purge_dirs; then
            rm -r output
            rm -r proc_data
        fi

        python proc_add_coord.py

    fi

    if $process_absem; then

        cd $date_dir/absem

        if $purge_dirs; then
            rm -r output
            rm -r proc_data
        fi

        python absem_setup.py

    fi


done