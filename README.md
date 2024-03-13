# MHD Photoionization

Modeling and experimental analysis of photoionization in MHD lab. 

## Installation 


a `.env` file needs to be created in the root repository directory with the following info
```
REPO_DIR='C:/Path/to/this/folder'
```

Initialize the mhdpy submodule

```
git submodule init
git submodule update
```

Also need to go in mhdpy and run these to the the 'readTrc' submodule for lecroy oscilloscope data. 

create a python virtual environment, using a python 3.10 installation

`python -m venv venv`

Then activate this envrionment and run `install.sh`


### Data Pipeline 

Notes on data pipeline. Processing from raw data on WSL. 

Data is on external hard drive D:/

1. create folder paths to H drive and Lecroy raw data in .env file (TODO: standardize names)
2. mount the D drive in wsl
    ```
    mkdir /mnt/d
    sudo mount -t drvfs D: /mnt/d
    ```