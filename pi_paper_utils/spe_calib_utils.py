
from skimage import transform
from pint import Quantity


# Scale determined from calbration image after applying the perspective transformation
CORRECTED_IMG_SCALE = Quantity(6/403, 'inch/pixel').to('mm/pixel').magnitude # 403 pixels per 6 inch

def transform_projective(ds, tform):

    # ds shoudl be a a dataset with x and y coordinates from 0 to 1024
    # da_out = ds_int['counts']

    ds_img = ds['counts'].values

    ds_img_tf = transform.warp(ds_img, tform.inverse)
    
    # need to flip the image vertically to line up correctly with coordinates
    ds['counts_corr'] = (('y', 'x'), ds_img_tf[::-1, :])
    # ds['counts_corr'] = (('y', 'x'), ds_img_tf)

    da_out = ds['counts_corr']

    return da_out

import numpy as np

def pipe_transform_projective(ds_sel, tform, downsel_range={'x': slice(0, 1024), 'y': slice(412, 612)}):

    # ds_sel is an image dataset with gatedelay coordinates and raw pixel coordinates 
    x_grid = np.arange(0, 1024, 1)
    y_grid = np.arange(0, 1024, 1)
    ds_sel_int = ds_sel.interp(x=x_grid, y=y_grid)#.fillna(0)

    # Squeeze added to remove the gatedelay dimension. new to 2024 xarray. https://github.com/pydata/xarray/pull/9280
    da_tf = ds_sel_int.groupby('gatedelay').apply(lambda x: transform_projective(x.squeeze(), tform))

    if downsel_range is not None:
        da_tf = da_tf.sel(downsel_range)

    da_tf.coords['x'] = (da_tf.coords['x'] - LOC_BARREL_EXIT) * CORRECTED_IMG_SCALE
    da_tf.coords['y'] = (da_tf.coords['y'] - LOC_BARREL_CENTERLINE) * CORRECTED_IMG_SCALE

    return da_tf


#Simple calibration
# Based on un perspective corrected image

CALIB_SIMPLE = Quantity(66.3, 'pixels/inch').to('pixels/mm').magnitude
LOC_BARREL_CENTERLINE = 512
LOC_BARREL_EXIT = 250

def calibrate_da_pimax_simple(da):
    #TODO: implement calbiration that takes into account perspective
    da['y'] = da['y'] - LOC_BARREL_CENTERLINE
    da['x'] = da['x'] - LOC_BARREL_EXIT

    da['x'] = da['x']/CALIB_SIMPLE
    da['y'] = da['y']/CALIB_SIMPLE

    return da