"""
testing perspective calibration for topcam 
https://medium.com/@flcamarao/image-processing-using-python-homography-matrix-a5da44f3a57b
"""

#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()

from skimage import transform
from skimage.io import imread, imshow
import cv2
import pickle


#TODO: incorporate into data pipeline
fp_in = pjoin(REPO_DIR, 'experiment','data','manual', 'PI Max Calibration', 'Log_PI_TopCam_0 2023 May 18 11_13_47_process.tif')

img = imread(fp_in)
imshow(img)

# %%
# Convert the image from BGR to RGB for proper visualization with Matplotlib
# img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

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

# save the transformation matrix object

fp_out = pjoin(DIR_DATA_OUT, 'tform_projective.pkl')
with open(fp_out, 'wb') as f:
    pickle.dump(tform, f)


