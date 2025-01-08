"""Simple workaround to read REPO_DIR from .env file and output to tex file that can be imported into paper to defined repo_dir command"""

import os

# Read environment variable
repo_dir = os.getenv('REPO_DIR')


with open('SI/repodir.tex', 'w') as f:
    f.write(r'\newcommand{\repodir}{' + repo_dir + '}')
