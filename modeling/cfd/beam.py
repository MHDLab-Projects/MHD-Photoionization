"""

Summary
-------

Work in progress


"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import tri_mesh

n_dim = 3
_source = np.array([0.1,0.0,-0.1])
_target = np.array([0.2,0.0,0.1])

def new_beam_dataset(s=[0.0,1.0],n_rays=1):
    coords = {"s":np.array(s),"i_ray":np.arange(n_rays,dtype="int"),"i_dir":["x","y","z"]}
    return xr.Dataset(coords)

def vec_mag(x):
    return np.sum(x**2, axis=-1)**0.5

def vec_norm(x):
    return x/vec_mag(x)

def coordinate_tensor(e_1, e_2, method="cross"):
    """array of coordinate system basis vectors given primary (axis) direction
    and reference vector for tangent

    Parameters
    ----------
    e_1 : float array[3]
        1st direction vector
    e_2 : float array[3]
        2nd direction vector

    Returns
    -------
    ehat : float array [3,3]
        tensor contain coordinate system basis vectors
        ehat[i,:] - basis vector of direction i

    """
    #
    # create a coordinate system centered at the beam source
    # oriented to the target
    #
    e_1 = np.array(e_1)
    e_2 = np.array(e_2)

    ehat = np.zeros((n_dim,n_dim))
    ehat[0] = vec_norm(e_1)
    ehat[1] = vec_norm(e_2)

    if method == "cross":
        ehat[2] = np.cross(ehat[0],ehat[1])
        ehat[1] = np.cross(ehat[2],ehat[0])
    else:
        ehat[1] = ehat[1] - ehat[0]
        ehat[2] = np.cross(ehat[0],ehat[1])

    for i in range(n_dim):
        ehat[i] = vec_norm(ehat[i])

    return ehat

class Beam:
    """
    Summary
    -------
    A multi-wavelength beam from a source to a target

    Attributes
    ----------
    source : float[n_dim]
        position of the source
    target : float[n_dim]
        position of the target
    diam : float
        diameter of the beam
    ds : xr dataset
        data about beam

    axis : float[n_dim]
        beam axis
    plane : float[2][n_dim]
        direction vectors in perpidular plane


    """

    def __init__(self, source=_source, target=_target, diam=1.0, e1=[0.0,1.0,0.0]):

        self.source = source
        self.target = target
        self.diam = diam
        self.e1 = vec_norm(np.array(e1))
        self.ds = new_beam_dataset()
        self._plane = np.zeros((2,n_dim))
        self._plane[0] = e1
        self._length = 1.0
        self._axis = np.zeros(n_dim)
        self.update_data()

    def update_data(self):
        self._length = vec_mag(self.target-self.source)
        self._axis[:] = vec_norm(self.target-self.source)
        shat = coordinate_tensor(self.axis, self.e1)
        self._plane[0,:] = shat[1]
        self._plane[1,:] = shat[2]
        ds = self.ds
        ds["source"] = "i_dir", self.source
        ds["target"] = "i_dir", self.target
        ds.attrs["diameter"] = self.diam
        ds.attrs["axis"] = self._axis
        ds.attrs["plane"] = self._plane

    @property
    def L(self):
        return self._length

    @property
    def axis(self):
        return self._axis

    @property
    def plane(self):
        return self._plane

    def adapt(self, x):
        pass

    def translate(self, x):
        pass

    def mesh(self, n_s=10, n_r=4, n_sides=6, outer_points="regular"):
        """

        Parameters
        ----------
        n_s : int, optional
            number of points in beam direction The default is 10.
        n_r : TYPE, optional
            number of radial points. The default is 4.
        n_sides : TYPE, optional
            number of sides of polygon. The default is 6.
        outer_points : TYPE, optional
            DESCRIPTION. The default is "regular".

        Returns
        -------
        None.

        """
        x = np.linspace(0, self.L, n_s)
        cyl = tri_mesh.new_poly_cylinder(x=x,r=self.diam/2.0,n_r=4,n_sides=6,outer_points="regular")
        shat = coordinate_tensor(self.axis, self._plane[0])
        ds = cyl["ds"]
        ds["pos_beam"] = ds["pos"]*1.0
        ds["pos"][:] = np.matmul(ds["pos_beam"].to_numpy(),shat) + self.source

        self.add_cart_coords(ds, "pos", "{}")
        self.add_cart_coords(ds, "pos_beam", "{}_beam")

        self.ds = ds
        self.update_data()


    @classmethod
    def add_cart_coords(cls, ds, pos="pos", fmt="{}"):
        for i, n in enumerate("xyz"):
            name = fmt.format(n)
            ds.coords[name] = ("s", "ray"), ds[pos].data[:,:,i]

def adapt_extruded(s, v, pos, f_v, rtol=1e-2, max_iter=10, do_plot=False):
    """adaption on an extruded unstructured mesh in the direction of extrusion

    Parameters
    ----------
    s[i_s] : float array
        distance coordinate in extruded direction
    pos[i_s,i_r,i_d] : float array
        initial set of coordinate vectors for function evaluation
    v[i_s,i_r,i_v] : float array
        field variables
    f_v(pos) : function
        function to calculate field variable given a position

    Notes
    -----
    i_s : extruded direction
    i_r : ray index in perpendicular plane
    i_v : variable index

    Returns
    -------
    dictionary on adapted emsh
        s : distance coordinate
        pos : position
        v : field variables @ positions

    Notes:
        the normalization needs to be calculated after each refinement
        if the initial sampling does not contain the maxima

    """
    w_small = 1e-6
    w_vsmall = 1e-12
    n_s, n_r, n_var = v.shape
    w = np.zeros(n_var)
    v1 = v.max(axis=(0,1))
    v0 = v.min(axis=(0,1))
    #print("max",v1)
    #print("min",v0)
    w = abs(v1 - v0) + (v1 + v0)*w_small + w_vsmall

    if do_plot:
        fig, ax = plt.subplots()
    for i_iter in range(max_iter):
        v1 = v.max(axis=(0,1))
        v0 = v.min(axis=(0,1))
        w = abs(v1 - v0) + (v1 + v0)*w_small + w_vsmall

        dv = (v[1:] - v[:-1])/w
        ii = np.where(np.max(abs(dv), axis=1) > rtol)
        i_new = np.unique(ii[0])
        #print(i_iter, len(i_new), np.max(abs(dv), axis=1)[:,0])
        #print(i_iter, len(i_new), np.max(abs(dv), axis=1)[:,1])
        if len(i_new) < 1: break

        s_new = (s[i_new] + s[i_new+1])*0.5
        pos_new = (pos[i_new] + pos[i_new+1])*0.5
        if do_plot:
            ax.plot(s_new, pos_new[:,0,0],'o',label="{}".format(i_iter))
        v_new = f_v(pos_new)
        v_new = np.transpose(v_new, (1,0,2))
        s = np.concatenate([s,s_new])
        pos = np.concatenate([pos,pos_new])
        v = np.concatenate([v,v_new])
        ii = np.argsort(s)
        v = v[ii,:,:].copy()
        s = s[ii].copy()
        pos = pos[ii,:,:].copy()

    return {"s":s,"pos":pos,"v":v}


def adapt_beam_dataset(ds, adapt_vars, field_func, rtol=1e-2, max_iter=10):
    """adaption on an extruded unstructured mesh in the direction of extrusion

    Parameters
    ----------
    ds : xr.Dataset
        beam dataset
    adapt_vars : list of strings
        variable names to adapt
    field_func : function
        f(pos) -

    Returns
    -------
    dictionary on adapted emsh
        s : distance coordinate
        pos : position
        v : field variables


    """
    pos = ds["pos"].copy().to_numpy()
    n_s, n_ray, n_dim = pos.shape
    n_var = 2
    v = np.zeros((n_s,n_ray,n_var))
    for i_var, name in enumerate(adapt_vars):
        v[:,:,i_var] = ds[name]
    s = ds["s"].copy().to_numpy() ##

    adapt_out = adapt_extruded(s, v, pos, field_func, rtol=rtol, max_iter=max_iter)

    coords = {"s":adapt_out["s"], "ray":ds.ray, "i_dir":ds.i_dir}
    ds_new = xr.Dataset(coords)
    for i, name in enumerate(adapt_vars):
        ds_new[name] = ("s","ray"), adapt_out["v"][:,:,i]

    ds_new["pos"] = ("s","ray","i_dir"), adapt_out["pos"]
    #for i, name in enumerate(ds.i_dir.to_numpy()):
    #    ds_new.coords[name] = ("s", "ray"),  ds_new["pos"].data[:,:,i]

    x0 = ds_new["pos"] #.data[0,:,:]
    x1 = ds_new["pos"] #.data[-1,:,:]
    s_new = ds_new.s

    #print(x0)
    # not need computed in add_cart_coodinates
    out = x0*(1-s_new) + x1*s_new

    #ds_new["pos_beam"] = ("s","ray","i_dir"), x0*(1-s_new) + x1*s_new
    ds_new["pos_beam"] = ds_new["pos"]*0.0
    ds_new["pos_beam"][:,:,0] = s_new
    ds_new["pos_beam"][:,:,1] = ds["pos_beam"][0,:,1]
    ds_new["pos_beam"][:,:,2] = ds["pos_beam"][0,:,2]

    Beam.add_cart_coords(ds_new, "pos", "{}")
    Beam.add_cart_coords(ds_new, "pos_beam", "{}_beam")


    for name in ds:
        if not (name in ds_new):
            if "s" not in ds[name].coords:
                ds_new[name] = ds[name]
            else:
                try:
                    print(name)
                    ds_new[name] = ds[name].interp({"s":ds_new["s"]})
                except Exception as e:
                    print(e)
                    print(name, "not added")

    return ds_new