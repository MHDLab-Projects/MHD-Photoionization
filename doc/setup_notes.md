# Setup Notes

Various notes getting the repo and data pipeline working on various computers. 

### Files to copy

files currently needing to copy. create some sort of automation sync for this

experiment/data/manual
final/figures/setup_images
modeling/cfd/simulation_data

### Installation Process WSL

Documenting process of setting up repository. Vscode running WSL on same laptop codes were developed on. 
palette: Ctrl + Shift + P
TODO: update/merge with above

1. clone MHD-Photoionization repo (checkout dev branch)
2. palette: create python environment (venv). Python 3.11.2.  use reqs/requirements.wsl
3. run source install.sh
4. setup .env file (rename .env_example to .env)
5. open terminal in `automation` and run `source munge.sh` (TODO: vscode task)
    TODO: Full munging untested. Copying from sharepoint for now. 
6. run processed munged data task in vscode


Now installing on 7400

#### WSL setup process

TODO: redo this on fresh install. Should be able to upgrade to bookworm. then install `git python inkscape`

1. Install wsl debian. install wsl, python, jupyter extensions in vscode and connect to wsl. open terminal
2. sudo apt update, sudo apt upgrade


3. sudo apt install git python (this installs python 3.9 as python3... not correct)

Failed attempt to install python3.11 through deadsnakes:
4. sudo apt install software-properties-common (https://phoenixnap.com/kb/add-apt-repository-command-not-found-ubuntu)
5. sudo apt install gpg
5. sudo add-apt-repository ppa:deadsnakes/ppa (https://ubuntuhandbook.org/index.php/2022/10/python-3-11-released-how-install-ubuntu/)
6. Ended up uninstalling deadsnakes as authentication errors (think it is only for Ubuntu )

Switch from bullseye to bookworm
7. upgrade from bullseye to bookworm. Change bulleye-> bookworm in /etc/apt/sources.list
8. sudo apt update, sudo apt upgrade
9. had to sudo umount -l /lib/modules/ because of a usrmerge error
10 sudo apt update/upgrade
10. reboot/ reload

#### repostiory setup

clone repository. 
follow install sequence above
Mount the hard drive below before running munging 

cant do everything because of memory errors- just manually copying proc_data and final dataset

had to install inkscape : sudo apt install inkscape

Installing latex: 

latex workshop extension 
sudo apt-get install texlive-full 


#### Missing figures

list of missing figues with data scattered around. Move to final data. For now going to copy files in their file strucutre. 
Series of scripts that need to run: 

experiment/notebook/2018-11-20/munge_spe.py
experiment/notebook/2018-11-20/analysis_beam_timing.py

experiment/analysis/topcam/spe_perspective.py
experiment/analysis/topcam/spe_analysis_2023-05-24.py

modeling/viability/dataset/cantera_gen.py
modeling/viability/dataset/calc_P_zero.py
modeling/viability/dataset/calc_gamma_new.py

modeling/viability/figure_panels/gen_figs.sh


--For SI--


experiment/analysis/auto/expt_overview.py

experiment/analysis/mws/resampling

experiment/analysis/topcam/spe_analysis_2023-05-18.py

### Data Pipeline 

#TODO: Had to remount when reloading vscode window. Automate/fix. 

Notes on data pipeline. Processing from raw data on WSL. 

Data is on external hard drive D:/

1. create folder paths to H drive and Lecroy raw data in .env file (TODO: standardize names)
2. mount the D drive in wsl
    ```
    mkdir /mnt/d
    sudo mount -t drvfs D: /mnt/d
    ```


