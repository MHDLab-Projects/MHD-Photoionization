# Run all python scripts and echo their name to the terminal

# for i in $(ls *.py); do
#     echo $i
#     python $i
# done

python gamma_boundary_enhance.py
python gamma_panels.py
python gas_species.py
python UV_absorption.py