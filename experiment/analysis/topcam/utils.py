
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

