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

# %%
# Convert the image from BGR to RGB for proper visualization with Matplotlib
# img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

# Define the source points
src = np.array([251, 478,
                716, 474,
                251, 549,
                715, 541,
]).reshape((4, 2))

# Plot the image
plt.imshow(img)

# Plot the points on the image with red color
plt.scatter(src[:, 0], src[:, 1], color='red', marker='o')

# %%
dst = np.array([251, 477,
                715, 477,
                251, 547,
                715, 547,
]).reshape((4, 2))
tform = transform.estimate_transform('projective', src, dst)

#%%

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

ds['counts'].mean('gatedelay').mean('estime').plot()

#%%
ds.coords['x']

#%%

xgrid = np.arange(0, 1024, 1)
ygrid = np.arange(0, 1024, 1)

xgrid, ygrid = np.meshgrid(xgrid, ygrid)

xgrid