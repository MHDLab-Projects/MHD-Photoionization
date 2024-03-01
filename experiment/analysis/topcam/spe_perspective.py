"""
testing perspective calibration for topcam 
https://medium.com/@flcamarao/image-processing-using-python-homography-matrix-a5da44f3a57b
"""

#%%
from mhdpy.analysis.standard_import import *

from skimage import transform
from skimage.io import imread, imshow
import cv2


#TODO: incorporate into data pipeline
fp_in = pjoin(REPO_DIR, 'experiment','data','manual', 'PI Max Calibration', 'Log_PI_TopCam_0 2023 May 18 11_13_47_process.tif')

img = imread(fp_in)
imshow(img)

#%%
img.shape

# %%
# Convert the image from BGR to RGB for proper visualization with Matplotlib
# img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

# # Define the source points
# src = np.array([251, 478,
#                 716, 475,
#                 251, 550,
#                 715, 540,
# ]).reshape((4, 2))


# Define the source points no y distorion (ruler may not have been held flat)
# Measure the width of the ruler in pixels. then calculate the positions based on a flat ruler
# Tweaked final value unti output image ruler had consistent inch spacing
# This will scale the image but not correct for y distortion
# Tweaked center value to get laser off image to line up. Shouldn't affect the scaling/projection
center_value = 509
dy_barrel = 71
dy_horns = 66

y_barrel_1 = center_value - dy_barrel/2
y_barrel_2 = center_value + dy_barrel/2

y_horns_1 = center_value - dy_horns/2
y_horns_2 = center_value + dy_horns/2

src = np.array([251, y_barrel_2,
                716, y_horns_2,
                251, y_barrel_1,
                715, y_horns_1,
]).reshape((4, 2))

# # Extreme distorion for testing 
# src = np.array([251, 476,
#                 716, 496,
#                 251, 545,
#                 715, 565,
# ]).reshape((4, 2))

# Plot the image
plt.imshow(img)

# Plot the points on the image with red color
plt.scatter(src[:, 0], src[:, 1], color='red', marker='o')

# %%
dst = np.array([251, 476,
                715, 476,
                251, 545,
                715, 545,
]).reshape((4, 2))
tform = transform.estimate_transform('projective', src, dst)

#%%

tf_img = transform.warp(img, tform.inverse)
fig, ax = plt.subplots()
ax.imshow(tf_img)


#%%

tform.inverse

#%%k

# save the transformation matrix object

import pickle

fp_out = pjoin(DIR_DATA_OUT, 'tform_projective.pkl')
with open(fp_out, 'wb') as f:
    pickle.dump(tform, f)



#%%

fp_out = pjoin(DIR_DATA_OUT, 'img_projective.tif')

# Save the image

tf_img_out = tf_img*255
# Convert the image data to 8-bit unsigned integers
tf_img_out = tf_img_out.astype(np.uint8)

# Convert the image from RGB to BGR for proper saving with OpenCV
# tf_img_out = cv2.cvtColor(tf_img_out, cv2.COLOR_RGB2BGR)

cv2.imwrite(fp_out, tf_img_out)

# %%

fp_test_spe = pjoin(REPO_DIR, 'experiment', 'data', 'munged', '2023-05-24', 'spe', 'PI_topcam_4.cdf')

ds = xr.load_dataset(fp_test_spe)

ds = ds.mean('estime') #Needed for large file size

#%%

ds['counts'].mean('gatedelay').plot()


#%%

da_test = ds.isel(gatedelay=0)['counts']

xs = np.linspace(0, 1024, 1024)
ys = np.linspace(0, 1024, 1024)

da_test_int = da_test.interp(x=xs, y=ys).fillna(0)


img_test = da_test_int.values


plt.axhline(50, color='r')

imshow(img_test[460:550, 100:600])
# imshow(img_test)

#%%

img_test_tf = transform.warp(img_test, tform.inverse)


# zoom in 
plt.axhline(50, color='r')

imshow(img_test_tf[460:550, 100:600])
# imshow(img_test_tf)

#%%
# interpolate to 2d grid and use skimage.transform.warp to apply the transformation

from experiment.analysis.topcam.calib_utils import pipe_transform_projective

da_tf = pipe_transform_projective(ds, tform)


#%%

# da_tf.plot(col='gatedelay')



# %%


# Calibraiton image

# Munged calibration image. 


fp_in = pjoin(REPO_DIR, 'experiment','data','munged', '2023-05-18', 'spe','calibration', 'Log_PI_TopCam_0 2023 May 18 11_13_47.cdf')

ds_cal = xr.load_dataset(fp_in)
ds_cal = ds_cal.isel(frame=0).isel(estime=0) #repetitive one exposure

ds_cal = ds_cal/ds_cal.max('x').max('y') # Normalize 


ds_cal['counts'].plot(robust=True)
# %%

from experiment.analysis.topcam.calib_utils import transform_projective, pipe_transform_projective

# Assign a dummy gatedelay, as pipe_transform_projective expects it...TODO: fix this
ds_cal = ds_cal.assign_coords(gatedelay=[0])

da_cal_tf = pipe_transform_projective(ds_cal, tform)

da_cal_tf = da_cal_tf.isel(gatedelay=0)

# %%


# fig, axes = plt

# da_cal_tf.sel(y=slice(-15,15)).plot(robust=True, figsize=(20, 5), cmap='gray')

# for i in range(8):
#     offset = Quantity(i, 'in').to('mm').magnitude

#     plt.axvline(offset, color='r', alpha=0.3, linestyle='--')


# plt.xlabel('x (mm)')
# plt.ylabel('y (mm)')


# plt.savefig(pjoin(DIR_FIG_OUT, 'calibration_image_projective.png'))

# %%
fig, ax = plt.subplots()
ax.imshow(img)

#%%

import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(10, 10))

gs = gridspec.GridSpec(2, 2, height_ratios=[4, 1], width_ratios=[1, 1])
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])
ax2 = plt.subplot(gs[2:])

# Display img on the upper left subplot
ax0.imshow(img, cmap='gray')
ax0.set_title('Raw Image')

# Display tf_img on the upper right subplot
ax1.imshow(tf_img, cmap='gray')
ax1.set_title('Projective Transformation')

# Display da_cal_tf on the bottom subplot
da_cal_tf.sel(y=slice(-15,0)).plot(ax=ax2, vmin=0.1, vmax=0.3, cmap='gray')

for i in range(8):
    offset = Quantity(i, 'in').to('mm').magnitude
    ax2.axvline(offset, color='r', alpha=0.3, linestyle='--')

ax2.set_xlabel('x (mm)')
ax2.set_ylabel('y (mm)')

ax2.set_title('Transformed and Scaled')

plt.tight_layout()
plt.savefig(pjoin(DIR_FIG_OUT, 'calibration_image_projective.png'))