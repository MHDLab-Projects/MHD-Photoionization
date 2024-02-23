"""
Analysis of laser profile taken from thorlabs camera
output cropped grayscale calibrated image from imagej

https://imagej.net/imaging/spatial-calibration
"""

#%%

from mhdpy.analysis.standard_import import *

from tifffile import TiffFile

fp_in = 'output/image_18_crop.tif'

# Get the XResolution and YResolution tags


with TiffFile(fp_in) as tif:
    images = tif.asarray()
    for page in tif.pages:
        for tag in page.tags.values():
            if tag.name == 'XResolution':

                xres = tag.value
                break

xres = xres[0] / xres[1]

scale = 1/xres # mm/pixel

scale

#%%
            
from PIL import Image
import numpy as np
import xarray as xr


from PIL import Image
import numpy as np
import xarray as xr

# Open the image file with PIL
img = Image.open(fp_in)

# Initialize a list to hold all frames
frames = []
while True:
    try:
        # Try to get the next frame.
        img.seek(img.tell()+1)
        frames.append(np.array(img))
    except EOFError:
        # If EOFError is raised, there's no more frames.
        break

# Convert the list of frames to a 3D numpy array
img_array = np.array(frames)

# Convert the 3D numpy array to an xarray Dataset
ds = xr.DataArray(img_array, dims=['frame', 'y', 'x']).to_dataset(name='image')

# Apply the calibration to x and y axes

ds['x'] = ds['x'] * scale
ds['y'] = ds['y'] * scale


da = ds['image']
#%%

da.mean('frame').plot.imshow(cmap='gray')

#%%

da_mean = da.mean('x').mean('y')

da_mean.plot()

#%%
# keep only frames where the laser is on

da_on = da.where(da_mean > 10)

da_on.mean('frame').plot.imshow(cmap='gray')

# %%

da_on_mean = da_on.mean('frame')

hori_line = da_on_mean.sel(y=8, method='nearest')
hori_line.plot(label='horizontal line')
vert_line = da_on_mean.sel(x=5.5, method='nearest')
vert_line.plot(label='vertical line')

# find the positions along the lines where the intensity is half of the maximum
hori_half_max = hori_line.where(hori_line > hori_line.max()/2, drop=True)
hori_half_max.plot(color='r')

vert_half_max = vert_line.where(vert_line > vert_line.max()/2, drop=True)
vert_half_max.plot(color='r')

plt.legend()
# 
# %%

vert_half_max

#%%

delta_vert = vert_half_max['y'].values[-1] - vert_half_max['y'].values[0]
delta_hori = hori_half_max['x'].values[-1] - hori_half_max['x'].values[0]

print(f'Vertical FWHM: {delta_vert:.2f} mm')
print(f'Horizontal FWHM: {delta_hori:.2f} mm')

area = delta_vert * delta_hori

print(f'Area: {area:.2f} mm^2')

#%%

from pint import Quantity

Quantity(area, 'mm^2').to('cm^2')

#%%

