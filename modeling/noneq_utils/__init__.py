import os

# This is copied from mhdpy
def chdir_if_nb_render():
    """
    This is a utility designed for a workflow where analysis scripts are
    rendered into a 'nb_render' subdirectory. In that case the working directory
    of the Jupyter notebook kernel will have to be changed to the parent
    directory of the scripts such that relative file paths are resolved
    correctly. This function will change into the parent directory if the
    current file is in a directory called 'nb_render' 

    To use this function add the following line at the beginning of the python
    script (that will be render into nb_render folder)

    `from noneq_utils import chdir_if_nb_render; chdir_if_nb_render()`

    """

    folder_name = os.path.split(os.getcwd())[1]
    if folder_name == 'nb_render':
        print("Current folder is called 'nb_render', changing working directory to parent")
        os.chdir('..')