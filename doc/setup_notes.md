# Setup Notes

Various notes getting the repo and data pipeline working on various computers. 

## Overview

Originally was pursuing WSL but found some issues. So ended up getting full pipeline working with git bash (same setup on office comp)
    * Seems to be more prone to memory errors. I think using git bash in native windows can take advantage of memory caching. 
    * timestamp issues with lecroy munging, see lecroy timestamp seciton. Switched to trigger time attribute. Think this should relove the problem on WSL too, however have switched to multiprocessing for lecroy data to speed up time, and pretty sure this causes memory crashes, but should try again. 

## General installation process

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

### Files to copy

files currently needing to copy. create some sort of automation sync for this

experiment/data/manual
final/figures/setup_images
modeling/cfd/simulation_data

### Git bash setup

Add inkscape to windows system path. This is method used before with LDES viability paper, see that readme. It says issues were encountered with inkscape python overriding venv, but I didn't encounter that. 

-- latex on windows setup --

latex workshop extension

Install pandoc from github releases msi installer. Install for all users. (worked without setting up latex. )

texlive:i https://www.tug.org/texlive/windows.html


https://strawberryperl.com/

### WSL setup 

clone repository. 
follow install sequence above
Mount the hard drive below before running munging 


#TODO: Had to remount when reloading vscode window. Automate/fix. 
Data is on external hard drive D:/

1. create folder paths to H drive and Lecroy raw data in .env file (TODO: standardize names)
2. mount the D drive in wsl
    ```
    mkdir /mnt/d
    sudo mount -t drvfs D: /mnt/d
    ```


cant do everything because of memory errors- just manually copying proc_data and final dataset

had to install inkscape : sudo apt install inkscape

Installing latex: 

latex workshop extension 
sudo apt-get install texlive-full 


#### Setting up WSL
Documenting WSL setup process on Lattitude 7400

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


# Lecroy timestamp 

The lecroy timestamps originally came from the mtime property of the trc file. 

I cannot remember why this was chosen exactly, but it seemed to work when processing data on office comp from the z drive. However, when moving to processing on backup of raw data on external drive with a laptop, this failed and the timestamps were unreliable. C time also appeared to be unreliable. 

I recently pulled the lecroy file attributes to be attributes of the output xr.dataset. There is a trigger time in there which seems to agree with the mtime from the office comp within ~100 ms when spot checking for the 516 case. 

I'm not sure why I didn't use trigger time before. Going to test how things are working with trigger time. First making a commit with the ts_option flag set to mtime to produce the old results on the office comp. 
