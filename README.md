# MHD Photoionization

Modeling and experimental analysis of photoionization in MHD lab. 

## Installation 

This repository is tested to work with 
* VS Code on windows. 
    * Need VS Code python extension installed
* Python 3.12 through a virtual environment
* Terminal using Git Bash (Git for windows)
    * Select `Terminal: Select Default Profile` in VS Code command palette to switch terminal to git bash.
    * Git Bash Side Notes:
        * As of 01-2025 there is currently a bug with the python extension that the directory will not show correctly in Git Bash. The directory is actually correct and you can run scripts fine, but you can also type `cd .` to fix the display. https://github.com/microsoft/vscode-python/issues/23382
        * It is also possible to use Windows Subsystem for Linux (or linux directly), but by default this can run into memory issues during data munging whereas Git Bash seems to do fine (presumably caching excess memory to disk). 
        * In general A linux shell is needed for various `.sh` automation scripts, but perhaps these could be converted to powershell scripts. 

### Procedure:

`"CP" refers to the VS Code Command Palette shortcut (Ctrl+Shift+P)`

* Clone the repository. Do `(CP) Git: Clone` and enter the `web-url.git` for this repository, which can be found at the top of this page: ![repo-link](/doc/media/web-url.png). Then, select the directory in which to clone this repository; for example, `C:/Users/USERNAME/code`.

* Copy and rename the `.env_example` to `.env`. The directories in `.env` need to be updated to match your local filepaths. 
    * Change `REPO_DIR` to the directory in which you cloned this repository (for example, `C:/Users/USERNAME/code/MHD-Photoionization`). 
    * Change `CFD_RESULTS_DIR` to the directory in which you copied the simulation data. For example, `C:/Users/USERNAME/data/Photoionization/Simulations`
    * TODO: Finalize and explain how extra input data will be downloaded and imported

* Create python virtual environment `(CP) Python: Create Environment`. Select your Python 3 interpreter. If no interpreters appear, click the refresh icon in the top right.
![doc image](/doc/media/interpereter_select.png)

* For dependencies to install, select `requirements_gitbash.txt`. After installing `(CP) Python: Select Interpreter` and select .venv

* Select a python script (e.g. `setup.py`) to activate the vscode python extension. Open a terminal `ctrl + ~` and ensure the python .venv is activated. 
    * type `which python` to make sure it points to the path in which you cloned the repository.

* run `source install.sh` in the terminal (or run steps inside manually)


* Add in extra input data: The repository requires some external files to fully reproduce the final results. How this is being handled is in flux but for now there is a 'Input Files' folder that contains two directories. The data are categorized based on whether they are currently reproduced by the repository from the original raw lab data. Go within these directories, and copy the contents into the base level of the directory. 
    * Folder 1: `Extra (Not Reproduced)`: This folder contains a few various data files that were manually created. 
        * these should be revisited and added to the repository pipelines where possible, but some can't (e.g. pictures)
        * There is some data that is placed in `modeling/em-sim`. This data uses the FDTD package which doesn't have a pip install currently. But this data can be reproduced by separately creating a conda environment (see `modeling\em-sim\environment.yml`)
    * Folder 2: `Processed (Reproduced)`. This folder contains datasets that are produced by the repository (with access to Raw Data). These datasets are included here to avoid the long processing time, file size, etc. In particular the `experiment/data/munged` folder is large and takes a while to reproduce. Data munging is commented out in the main data pipeline (`automation/pipe_main.sh`) by default. 

Now the repository should be setup. 

There are various VS Code tasks that can be used to run the data processing (defined in `.vscode/tasks.json`,  `(CP) Tasks: Run Task -> continue without scanning output`). These mainly run shell scripts in the `automation` directory, which can also be done manually with e.g. `source script_name.sh`

### Generate Main Paper Data/Figures


The main paper data pipeline starts from the munged dataset and generates the final paper figures. 

Run the `Main Paper Pipeline` task

Final figures are output in `final\figures\output`

### Supplementary data

for supplementary information and notebooks run the `Supplementary Pipeline` task. 

Then render doc\SI\SI.tex with Texlive in VSCode
    Debug errors (e.g. missing files) by searching for error in  `doc\SI\SI.log`

