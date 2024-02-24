"""
Analysis of laser profile taken from thorlabs camera
output cropped grayscale calibrated image from imagej

https://imagej.net/imaging/spatial-calibration
"""

#%%

from mhdpy.analysis.standard_import import *

from tifffile import TiffFile

fp_in = pjoin(REPO_DIR, 'experiment','data','manual', 'Thorcam', 'image_18_crop.tif')

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

fig, axes = plt.subplots(1,2, figsize=(10,4))



plt.sca(axes[0])

vertical_xs = 5.5
horizontal_xs = 8
# keep only frames where the laser is on

da_on = da.where(da_mean > 10)

da_on.mean('frame').plot.imshow(cmap='gray')

plt.axhline(horizontal_xs, color='r', linestyle='--')
plt.axvline(vertical_xs, color='r', linestyle='-.')

axes[0].set_title('Image of laser beam')
axes[0].set_xlabel('x (mm)')
axes[0].set_ylabel('y (mm)')

plt.sca(axes[1])
da_on_mean = da_on.mean('frame')

hori_line = da_on_mean.sel(y=horizontal_xs, method='nearest')
hori_line.plot(label='horizontal line', linestyle='--')
vert_line = da_on_mean.sel(x=vertical_xs, method='nearest')
vert_line.plot(label='vertical line', linestyle='-.')

# find the positions along the lines where the intensity is half of the maximum
hori_half_max = hori_line.where(hori_line > hori_line.max()/2, drop=True)
hori_half_max.plot(color='r', alpha=0.5)

vert_half_max = vert_line.where(vert_line > vert_line.max()/2, drop=True)
vert_half_max.plot(color='r', alpha=0.5)

delta_vert = vert_half_max['y'].values[-1] - vert_half_max['y'].values[0]
delta_hori = hori_half_max['x'].values[-1] - hori_half_max['x'].values[0]

area = delta_vert * delta_hori

# Add the 'delta_vert', 'delta_hori', and 'area' values to the plot
plt.text(0.05, 0.95, f'Vertical FWHM: {delta_vert:.2f} mm', transform=plt.gca().transAxes)
plt.text(0.05, 0.90, f'Horiz. FWHM: {delta_hori:.2f} mm', transform=plt.gca().transAxes)
plt.text(0.05, 0.85, f'Area: {area:.2f} mm^2', transform=plt.gca().transAxes)

plt.ylim(0, 100)

axes[1].set_title('')
axes[1].set_xlabel('Position along line (mm)')

plt.legend()

plt.savefig(pjoin(DIR_FIG_OUT, 'laser_profile.png'), bbox_inches='tight')
# 
# %%
#%%

from pint import Quantity

Quantity(area, 'mm^2').to('cm^2')

#%%

