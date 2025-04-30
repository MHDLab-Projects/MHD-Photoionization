#
!/bin/bash
``
source ../../.env

nb_dir=$EXPT_NB_DIR

dates=("2023-04-07" "2023-05-12" "2023-05-18" "2023-05-24")

process_mwt=true
process_aas=true

purge_dirs=false

for date in ${dates[*]};

do
    echo --------$date---------

    date_dir=$nb_dir/$date

    if $process_mwt; then

        cd $date_dir/mwt

        if $purge_dirs; then
            rm -r output
            rm -r proc_data
        fi

        python proc_add_coord.py

    fi

    if $process_aas; then

        cd $date_dir/aas

        if $purge_dirs; then
            rm -r output
            rm -r proc_data
        fi

        python aas_setup.py

    fi


done