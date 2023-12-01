date_list=('2023-04-07' '2023-05-12' '2023-05-18' '2023-05-24')

for date in "${date_list[@]}"
do
    python absem_setup.py -d $date
    # python combine_lecroy_time.py -d $date
done

python collect_data.py
python proc_add_coord.py