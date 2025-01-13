from mhdlab.analysis.standard_import import *
create_standard_folders()

from skimage import transform
from skimage.io import imread, imshow
import cv2
import pickle
from pi_paper_utils.spe_calib_utils import pipe_transform_projective

DATA_DIR_OUT = pjoin(REPO_DIR, 'final','dataset','topcam','output')
fp_transform = pjoin(DATA_DIR_OUT, 'tform_projective.pkl')
with open(fp_transform, 'rb') as f:
    tform = pickle.load(f)

#%%



#TODO: incorporate into data pipeline
fp_in = pjoin(REPO_DIR, 'experiment','data','manual', 'PI Max Calibration', 'Log_PI_TopCam_0 2023 May 18 11_13_47_process.tif')

img = imread(fp_in)
imshow(img)


tf_img = transform.warp(img, tform.inverse)
fig, ax = plt.subplots()
ax.imshow(tf_img)

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