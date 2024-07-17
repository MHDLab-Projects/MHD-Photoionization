# MHD Photoionization

Modeling and experimental analysis of photoionization in MHD lab. 

## Installation 

Rename `.env_example` -> `.env`. Update values as needed.  #TODO: remove unused variables/clean up. 

create python virtual environment using one the requirement files in 'reqs'

Tested with Python 3.10.11 and git bash. 

run `source install.sh` or run steps inside

copy in extra input data (folder of manual dataset files)
    must happen after submodule cloning (for em-sim data files. )

`cd automation`
`source full_pipeline.sh`


Then render doc\SI_man\SI_man.tex with Texlive in VSCode
    Debug errors (e.g. missing files) by searching for error in  `doc\SI_man\SI_man.log`

