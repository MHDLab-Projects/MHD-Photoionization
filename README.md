# MHD Photoionization

Modeling and experimental analysis of photoionization in MHD lab. 

## Installation 

Rename `.env_example` -> `.env`. Update values as needed.  

create python virtual environment using one the requirement files in 'reqs'

Tested with Python 3.12.1 and git bash. 

run `source install.sh` or run steps inside

copy in extra input data (folder of manual dataset files)
    EM simulation in 'modeling\em-sim' is a submodule. This uses meep FDTD package which does not have a pip install and simulations are run separately in a conda environment. 
    must happen after submodule cloning
    

`cd automation`
`source pipe_main.sh`

for supplementary information and notebooks

`source pipe_supplementary.sh`

Then render doc\SI_man\SI_man.tex with Texlive in VSCode
    Debug errors (e.g. missing files) by searching for error in  `doc\SI_man\SI_man.log`

