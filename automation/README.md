Various automation shell scripts

These scripts are designed to be run as commands in vscode tasks (defined in .vscode/tasks.json)

# Setup

need to give execute permission for vscode to run  

`chmod +x ./automation/render_all_expt.sh`



# vscode tasks.json

vscode tasks seems to run the scripts from the repository directory. therefore the .env file is just loaded with `source .env` instead of `source ../.env` but to run the scripts stanalone the latter could be used. 