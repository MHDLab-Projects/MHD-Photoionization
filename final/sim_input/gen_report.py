#%%
from mhdlab.analysis.standard_import import *
import re

output_dir = pjoin('output')

fname = pjoin(output_dir, 'sim_input_report.tex')
#%%


import collections

signal_stat_image_dir = pjoin(DIR_FIG_OUT, 'signal_stats_repeat')
# signal_stat_image_dir = os.path.abspath(signal_stat_image_dir)

# Initialize an empty dictionary using defaultdict
signal_images = collections.defaultdict(list)

for image_fn in os.listdir(signal_stat_image_dir):
    # filenames have format 'kwt_signalstr.png, extract kwt and signalstr. signalstr has underscores 
    regex = r'(\w+?)_(.+).png'

    matches = re.findall(regex, image_fn)

    kwt = matches[0][0]
    signal_str = matches[0][1]

    # Append the image filename to the list of images for this signal_str
    signal_images[signal_str].append(image_fn)

#%%
import collections
import os
import re

individual_signals_dir = os.path.join(DIR_FIG_OUT, 'individual_signals')

# Initialize an empty dictionary using defaultdict
signal_time_images = collections.defaultdict(list)

for image_fn in os.listdir(individual_signals_dir):
    # filenames have format 'date_signalstr.png', extract date and signalstr. signalstr has underscores 
    regex = r'(\d+?)_(.+).png'

    matches = re.findall(regex, image_fn)

    if matches:
        date = matches[0][0]
        signal_str = matches[0][1]

        # Append the image filename to the list of images for this signal_str
        signal_time_images[signal_str].append(image_fn)


#%%

from PIL import Image as PILImage
def get_figure_aspect(image_fp):
    pil_image = PILImage.open(image_fp)
    width, height = pil_image.size
    aspect = height / float(width)
    return aspect

#%%
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, PageTemplate, NextPageTemplate
from  reportlab.platypus.tableofcontents import TableOfContents
from  reportlab.platypus.frames import Frame
from reportlab.lib.styles import getSampleStyleSheet

pdf_output_path = pjoin('output', 'sim_input_report.pdf')

doc = SimpleDocTemplate(pdf_output_path, pagesize=letter)
story = []

styles = getSampleStyleSheet()

def _create_title_page(canvas, doc):
    canvas.saveState()
    canvas.setFont('Times-Bold', 16)
    canvas.drawCentredString(letter[0] / 2, letter[1] - 50, "Simulation Input Report")
    canvas.restoreState()

def _create_toc_page(canvas, doc):
    canvas.saveState()
    canvas.setFont('Times-Bold', 16)
    canvas.drawCentredString(letter[0] / 2, letter[1] - 50, "Table of Contents")
    canvas.restoreState()

# Create a different frame for the first page (title page)
frame1 = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
template1 = PageTemplate(id='Title', frames=frame1, onPage=_create_title_page)
doc.addPageTemplates([template1])

# Create a different frame for the second page (table of contents)
frame2 = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
template2 = PageTemplate(id='TOC', frames=frame2, onPage=_create_toc_page)
doc.addPageTemplates([template2])

# Create a normal frame for the rest of the pages
frame3 = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
template3 = PageTemplate(id='normal', frames=frame3)
doc.addPageTemplates([template3])


# Add a title page
story.append(NextPageTemplate('Title'))
story.append(Paragraph('Simulation Input Report', styles['Title']))
story.append(Paragraph('This report contains plots of the experimental signals used to generate the CFD simulation input file', styles['Normal']))
story.append(PageBreak())


#TODO: linked table of contents
# # Add a table of contents
# toc = TableOfContents()
# story.append(NextPageTemplate('TOC'))
# story.append(toc)
# story.append(PageBreak())

story.append(NextPageTemplate('normal'))

image_fp = pjoin(DIR_FIG_OUT, 'sim_input_timewindows.png')

aspect = get_figure_aspect(image_fp)

img = Image(image_fp, width=300, height=(300 * aspect))  # Maintain aspect ratio

story.append(img)

caption = 'Test case time windows used for simulation report. green= test case time windows'
story.append(Paragraph('<strong>Figure:</strong> %s' % caption, styles['Normal']))

story.append(PageBreak())


# add signal time images

for signal in signal_time_images:
    image_fns = sorted(signal_time_images[signal])
    story.append(Paragraph('Signal Time Series for signal: {}'.format(signal), styles['Heading1']))
    for image_fn in image_fns:
        date = image_fn.split('_')[0]
        img_path = pjoin(individual_signals_dir, image_fn)

        aspect = get_figure_aspect(img_path)

        img = Image(img_path, width=300, height=(300 * aspect))  # Maintain aspect ratio
        story.append(img)

        caption = 'signal: {}, date {}'.format(signal, date)
        story.append(Paragraph('<strong>Figure:</strong> %s' % caption, styles['Normal']))

        story.append(Spacer(1, 12))

    story.append(PageBreak())


# add signal statistic images
for signal in signal_images:
    image_fns = sorted(signal_images[signal])
    story.append(Paragraph('Signal Statistics for signal: {}'.format(signal), styles['Heading1']))
    for image_fn in image_fns:
        kwt = image_fn.split('_')[0]
        img_path = pjoin(signal_stat_image_dir, image_fn)

        aspect = get_figure_aspect(img_path)

        img = Image(img_path, width=300, height=(300 * aspect))  # Maintain aspect ratio
        # img.caption = 'kwt: {}'.format(kwt.replace('p', '.'))

        story.append(img)
        caption = 'signal: {},  kwt {}'.format(signal, kwt.replace('p', '.'))
        story.append(Paragraph('<strong>Figure:</strong> %s' % caption, styles['Normal']))
        story.append(Spacer(1, 12))
    story.append(PageBreak())

doc.build(story)