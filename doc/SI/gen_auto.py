import numpy as np

from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Matrix, Alignat, NoEscape, Command

from pylatex.utils import italic
import os
from os.path import join as pjoin

import dotenv; dotenv.load_dotenv()

if __name__ == '__main__':

    geometry_options = {"tmargin": "1cm", "lmargin": "2cm"}
    doc = Document(geometry_options=geometry_options)

    fig_folder = pjoin(os.getenv('REPO_DIR'), 'experiment','analysis','auto','output', 'absem_1d', '53x')

    doc.preamble.append(Command('title', 'Auto Generated Figures'))
    doc.append(NoEscape(r'\maketitle'))


    doc.append("""This document was generated automatically from the absem1d folder in the auto folder of the experiment/analysis directory. It contains all the figures in that folder.""")

    doc.append(NoEscape(r'\tableofcontents'))

    with doc.create(Subsection('Absorption Emission')):
        doc.append('Absorption Emission Figures')

        for fn in os.listdir(fig_folder):
            if fn.endswith('.png'):
                # doc.append(fn)
                image_filename = pjoin(fig_folder, fn)
                with doc.create(Figure(position='p')) as fig:
                    fig.add_image(image_filename, width=NoEscape(r'0.8\textwidth'))
                    fig.add_caption(fn)

        doc.append(NoEscape('\clearpage')) # clearpage forces output of stored figures and tables. 

    
    with doc.create(Subsection("New section")):
        doc.append('Some text.')

    doc.generate_pdf('auto', clean_tex=True)
    # doc.generate_tex('auto')