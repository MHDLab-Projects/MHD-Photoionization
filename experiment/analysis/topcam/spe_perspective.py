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
center_value = 512
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

ds['counts'].mean('gatedelay').mean('estime').plot()



#%%

ds_sel = ds.mean('estime')


#%%

#%%

da_test = ds_sel.isel(gatedelay=0)['counts']

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

da_tf = pipe_transform_projective(ds_sel, tform)


#%%

da_tf.plot(col='gatedelay')




# %%

# plt.xlim(400, 600)



# %%
#%%

# Trying to apply the transformation to the xarray dataset


# ds_stack = ds_sel.stack(temp=['x','y'])

# ds_stack

# mat = tform.inverse.params

# def transform_projective(ds, mat):
#     x = ds['x'].item()
#     y = ds['y'].item()
#     ds = ds.drop('temp')

#     x_new, y_new, z = np.dot(mat, [x, y, 1])

#     ds['x'] = x_new
#     ds['y'] = y_new

#     ds = ds.expand_dims('x').expand_dims('y').stack(temp=['x','y'])

#     return ds


# ds_stack.groupby('temp').apply(lambda x: transform_projective(x, mat))
#