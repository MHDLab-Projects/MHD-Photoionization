import numpy as np
import xarray as xr

def xr_add_pv(ds, pv_data, verbose=0):
    """add pyvista arrays from *data* to xarray Dataset, *ds*
    
    Parameters
    ----------
    ds : xarray Dataseries
        
    pv_data : pyVista array
    
    verbose : int, default 0
        how much data to print to screen
    
    
    """
    coord_names = [x for x in ds.coords]
    
    dim = coord_names[0]
    dirs = ds[coord_names[1]].values
    # print("dim = ", dim, "dir", dirs)
    for name, v in pv_data.items():
        if verbose > 0: 
            print(name, v.shape)
        elif len(v.shape) == 1:
            ds[name] = dim, v
        elif len(v.shape) == 2:
            for i, n in enumerate(dirs):
                ds[name+"_"+n] = dim, v[:,i]
        elif len(v.shape) == 3:
            for i, n in enumerate(dirs):
                for j, m in enumerate(dirs):
                    ds[name+"_"+n+m] = dim, v[:,i,j]
        else:
            print("rank-3 tensors not supported")
            
 
            
def pv_to_xr(mesh, flat=True):
    """create xarrays from pvista unstructured mesh data
    
    Parameters
    ----------
    mesh : pyVista mesh
    
    flat : bool, default True
        true  : store coordinates as separate 1D arrays
        false : store coordinates as a 2D array eg. (n_points, n_dim)
    
    Returns
    -------
    dict
        point : xarray for cell verticies
        cell  : xarray for cell averages - position at cell centers
        
    Note:
        flat is only implemented for the coordinates not other vectors (e.g. velocity, electric field, ... )
    
    """
    
    def add_points(ds, pos):
        coord_names = [x for x in ds.coords]        
        dim = coord_names[0]
        if flat:
            for i, n in enumerate("xyz"):
                ds["pos_"+n] = dim, pos[:,i] 
        else:
            ds["pos"] = (dim,"dir"), pos
    
    m = mesh
    ds_point = xr.Dataset(coords={"i_point":np.arange(m.n_points), "dir":["x","y","z"]})
    add_points(ds_point, m.points)
    xr_add_pv(ds_point, m.point_data)
    
    
    ds_cell =  xr.Dataset(coords={"i_cell":np.arange(m.n_cells), 
                                  "dir1":["x","y","z"], "dir2":["x","y","z"] })
    add_points(ds_cell, m.cell_centers().points )    
    xr_add_pv(ds_cell, m.cell_data)
   
    
    return {"point": ds_point, "cell":ds_cell}

import pyvista as pv

#TODO: No straightforward way to do this?
def downsel_arrays(mesh, keys):
    for k in mesh.cell_data.keys():
        if k not in keys: mesh.cell_data.remove(k)

    for k in mesh.point_data.keys():
        if k not in keys: mesh.point_data.remove(k)

    return mesh

def pv_to_unstack_xr(mesh):

    ds = pv_to_xr(mesh)['point']
    ds = ds.drop('dir')

    ds2 = ds.set_index(i_point=('pos_x', 'pos_y','pos_z'), append=True)
    ds2 = ds2.reset_index('i_point_level_0', drop=True).unstack('i_point')

    # ds2 = ds2.rename({
    #     'pos_x': 'x',
    #     'pos_y': 'y',
    #     'pos_z': 'z',
    # })

    return ds2