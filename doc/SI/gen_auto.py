import numpy as np

from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Matrix, Alignat, NoEscape
from pylatex.utils import italic
import os
from os.path import join as pjoin

import dotenv; dotenv.load_dotenv()

if __name__ == '__main__':

    geometry_options = {"tmargin": "1cm", "lmargin": "2cm"}
    doc = Document(geometry_options=geometry_options)

    fig_folder = pjoin(os.getenv('REPO_DIR'), 'experiment','analysis','auto','output', 'absem_1d', '53x')

    with doc.create(Subsection('Absem Auto')):
        doc.append('Absorption Emission')

        for fn in os.listdir(fig_folder):
            if fn.endswith('.png'):
                doc.append(fn)
                image_filename = pjoin(fig_folder, fn)
                with doc.create(Figure(position='p')) as fig:
                    fig.add_image(image_filename, width=NoEscape(r'0.8\textwidth'))
                    fig.add_caption(fn)

    doc.generate_pdf('auto', clean_tex=False)
    # doc.generate_tex('auto')