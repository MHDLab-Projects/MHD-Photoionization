"""
Summary
-------
Functions to do beam integration using xarray datasets.


See also xrBeamExample.ipynb

"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import pyvista as pv
import triangle

from pv_axi_utils import pv_interpolator, AxiInterpolator
import tri_mesh
from scipy import interpolate

n_dim = 3

def test_intensity_interpolator(ds, n_wave=30, n_extra=20):
    """Beam intensity interpolation interpoation testing

    Notes
    -----
    Idea - use a "wide-band" method for wavelength discretization

    1. take experimental intensity vs. wavelength function as create
        an interpolating polynomial of the function

    2. perform beam integration over wavelength interval to reduce the number
        of wavelength which need to be integrated

    work in progress
    """
    a = ds["diff"].sel(mp="barrel").isel(time=0)

    wave = np.linspace(a.wavelength.min(),a.wavelength.max(), 1000)
    data = a.data
    dwave = (wave[1] - wave[0])*n_extra

    n = n_wave
    wave_in = np.concatenate( [ np.linspace(wave[0]-dwave,wave[0],n+1)[:-1], a.wavelength, np.linspace(wave[-1], wave[-1]+dwave,n+1)[1:]] )
    data_in = np.concatenate( [ data[0]*np.ones(n), data, data[-1]*np.ones(n)]  )

    fig, axArray = plt.subplots(2,2)
    ax = list(axArray.flat)
    for ax_i in ax[:-1]:
        ax_i.plot( wave_in, data_in, '.', alpha=0.1)
    max_deriv = []
    for s in [None,0,1,10,100,1e3,1e4]:
        f_source = interpolate.UnivariateSpline(wave_in, data_in, ext=3, s=s)
        label = s
        if s == None:
            label = "none"
        for i, ax_i in enumerate(ax[:-1]):
            ax_i.plot(wave,f_source(wave),'-', alpha=0.5, label=label)

    ax[-1].plot(wave, f_source.derivative(1)(wave))


def laplacian_filter(x, u_in, alpha=0.5, n_filt=10):
    """1D laplacian filter

    Parameters
    ----------
    x : float array []
        grid spacing
    u : float array []
        array to filter
    alpha : float < 1
        strength of the filter


    """
    u = u_in.copy()
    dx = np.diff(x)
    dx_c = x*1.0
    dx_c[1:-1] = 0.5*(dx[1:] + dx[:-1])
    dx_c[0] = dx[0]
    dx_c[-1] = dx[-1]
    D = 0.5*alpha*dx/(1.0/dx_c[1:] + 1.0/dx_c[:-1])
    D = D.min()
    for i in range(n_filt):
        du = D*np.diff(u)/dx
        u[1:-1] += 2.0*np.diff(du)/(dx[1:] + dx[:-1])
    return u


def filter_wave_profile(ds, alpha=0.5, n_filt=1000):
    """apply 1D laplacian filter to ds

    Parameters
    ----------
    x : float array []
        grid spacing
    u_in : float array []
        array to filter
    alpha : float < 1
        strength of the filter
    """
    ds_filt = ds * 1.0
    x = ds_filt.wavelength
    mp_list = ds.mp
    time_list = ds.time


    for j, t in enumerate(time_list):
        for i, mp in enumerate(mp_list):
            u_all = ds_filt.sel(time=t, mp=mp)
            for name, u in u_all.items():
                v = laplacian_filter(x.to_numpy(), u.to_numpy(), alpha=alpha, n_filt=n_filt)
                ds_filt[name].loc[{"time":t, "mp":mp}]= v

    return ds_filt

#
# read wavelength dependence of source intensity
#
data_path = r"C:\Users\Huckaby\Desktop\MHD\Simulations\viz_example"
fname = "ds_calib.cdf"
source_intensity_fname = os.path.join(data_path)
def read_intensity_wavelength_function(fname=source_intensity_fname, do_plot=True):
    ds = xr.open_dataset(fname)
    n_vars = len(ds)
    n_pos = len(ds.mp)

    if do_plot:
        fig, ax = plt.subplots(n_vars,n_pos,sharex=True)

        for j, mp in enumerate(np.array(ds.mp)):
            ds1 = ds.sel(mp=mp).isel(time=0)
            for i, var_name in enumerate(ds1):
                ax[i,j].plot(ds1[var_name]["wavelength"], ds1[var_name], label="{}".format(var_name))
                ax[0,j].set_title(mp)
                ax[-1,j].set_xlabel("wavelength")
                ax[i,j].legend(loc='upper right')
        fig.tight_layout()
    return ds

#
# read solution and create interpolator
#
viz_path = r"C:\Users\Huckaby\Desktop\MHD\Simulations\viz_example"
fname = "frontCyl_plasma.vtk"
soln_fname = os.path.join(viz_path, fname)
def new_soln_interpolator(fname=soln_fname, do_plot=True):
    reader = pv.get_reader(fname)
    soln = reader.read()
    f_axi = AxiInterpolator(soln, var_names="T K".split())


    x_exit = 207.910 * 1e-3
    if do_plot:

        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        fig, ax = plt.subplots()
        center_line = soln.sample_over_line([x_exit,0,0],[0.3,0,0])
        ax.plot(center_line.points[:,0], center_line["T"])
        ax.set_ylabel("T [K]")
        tw = ax.twinx()
        for i, var in enumerate("K KOH".split()):
            tw.plot(center_line.points[:,0], 100*center_line[var],'--',label=var, color=colors[i+1])
        tw.legend()
        tw.set_ylabel("wt%")
        x_center = x_exit + np.array([20.0,25.0,30.0,35,40.0,45.0,50.0])*1e-3
        pos_center = np.zeros((len(x_center),3))
        pos_center[:,0] = x_center
        tw.plot(x_center, 100.0*f_axi(pos_center)[:,1],'o')

    return f_axi

def adapt_xr(ds_in, f_val, max_iter, rtol=0.01, verbose=0):

    ds_out = ds_in.copy()

    for i in range(max_iter):
        dv_ds = ds_out.diff("s")

def adapt_array(s, v, pos, f_v, max_iter=10, rtol=1e-2, do_plot=False):
    """one-dimensional adaption on an extruded unstructured mesh

    Parameters
    ----------
    s[i_s] : float array
        distance coordinate
    pos[i_s,i_r,i_d] : float array
        coordinate vector - for function evaluation
    v[i_s,i_r,i_v] : float array
        field variable
    f_v(pos) : function
        output variable vector

    Returns
    -------
    dictionary on adapted emsh
        s : distance coordinate
        pos : position
        v : field variables

    Notes:
        may not work if s, pos, v are xr.Dataarray, since the index values
        have a coordinate associated


    TODO: Need to

    """
    w_small = 1e-6
    w_vsmall = 1e-12
    n_s, n_r, n_var = v.shape
    w = np.zeros(n_var)
    v1 = v.max(axis=(0,1))
    v0 = v.min(axis=(0,1))
    w = (v1 - v0) + (v1 + v0)*w_small + w_vsmall

    if do_plot:
        fig, ax = plt.subplots()
    for i_iter in range(max_iter):
        dv = (v[1:] - v[:-1])/w
        ii = np.where(np.max(abs(dv), axis=1) > rtol)
        i_new = np.unique(ii[0])
        if len(i_new) < 1: break
        s_new = (s[i_new] + s[i_new+1])*0.5
        sA = s[i_new]*0.5
        sB = s[i_new+1]*0.5
        sC = sA[:] + sB[:]
        pos_new = (pos[i_new,:,:] + pos[i_new+1,:,:])*0.5
        pA = pos[i_new,:,:]
        pB = pos[i_new+1,:,:]
        pC = (pA + pB)/2.0
        if do_plot:
            ax.plot(s_new[:], pos_new[:,0,0],'o',label="{}".format(i_iter))
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



def adapt_beam(ds, f_axi, adapt_vars=["T","K"], do_plot=True):
    """

    Parameters
    ----------
    ds : xr.Dataset
        beam dataset
    f_axi : TYPE
        adaptation function
    adapt_vars : list of floats, optional
        variables to use for adapation functions, The default is ["T","K"].

    Returns
    -------
    ds_new : xr.Dataset
        new dataset after refinement

    """
    pos = ds["pos"].copy().to_numpy()
    n_s, n_ray, n_dim = pos.shape
    n_var = 2
    v = np.zeros((n_s,n_ray,n_var))
    for i_var, name in enumerate(adapt_vars):
        v[:,:,i_var] = ds[name]
    s = ds["s"].copy().to_numpy()
    out = adapt_array(s, v, pos, f_axi, do_plot=do_plot)

    coords = {"s":out["s"], "ray":ds.ray, "i_dir":ds.i_dir}
    ds_new = xr.Dataset(coords)
    for i, name in enumerate(adapt_vars):
        ds_new[name] = ("s","ray"), out["v"][:,:,i]

    for i, name in enumerate(ds.i_dir.to_numpy()):
        ds_new.coords[name] = ("s", "ray"), out["pos"][:,:,i]

    ds_new["pos"] = ("s","ray","dir"), out["pos"]

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


def calc_absorption(source, beam_ds, Q_abs):
    """calculate and add absorption to dataset

    Parameters
    ----------
    source : xr.Dataset
        contains total intensity at source [W]
    beam_ds : xr.Dataset
        beam Dataset
    Q_abs : list of Gaussian Absorption Cross-sections

    adds the following values to the beam
        I_source : intensity at the source [W/m^2]
        Q        : spectral absorption at each point [m^2]
        kappa    : absorption coefficient [1/m]
        kappa_kL : integrated abs. coefficient []
        I_target : intensity at the target [W/m^2]
        I_target_sum : total intensity at the target [W]

    """
    ds = beam_ds

    lam_nm = source.wavelength
    t = (ds["T"]/Q_abs.T0)
    #print(t.shape)

    #
    #  evaluate cross-section - source be moved to absorption function
    #   Q_abs.calc(T)
    #
    beta_list = [beta*t**Q_abs.b[i] for i, beta in enumerate(Q_abs.beta)]
    #print(beta_list[0])
    #val_list = [ beta*(lam_nm-Q_abs.lam[i])**2 for i,beta in enumerate(beta_list)]

    # cross-section of both peaks
    Q_list = [ Q_abs.A[i]*np.exp( -(beta*(lam_nm-Q_abs.lam[i]))**2) for i,beta in enumerate(beta_list)]

    ds.coords["wavelength"] = lam_nm
    # spectral obsorption as each point
    ds["Q"] = sum(Q_list)

    # intensity at the source
    ds["I_source"] = source["diff"].isel(time=0).drop("time")

    # absorption coefficient
    ds["kappa"] = ds["Q"]*ds["K"]

    ii = np.arange(len(ds.s))
    s = ds["s"].to_numpy()
    kappa = ds["kappa"].to_numpy()
    kappaL = kappa*0.0
    kappaL[1:] = (kappa[1:] + kappa[:-1])*(s[1:] + s[:-1])[:,np.newaxis,np.newaxis]*0.5
    ds["kappa_dL"] = ("s","ray","wavelength"), kappaL
    ds["kappa_L"] = ds["kappa_dL"].sum(dim="s")
    ds["I_target"] = ds["I_source"]*np.exp(-ds["kappa_L"])
    ds["I_target_sum"] = ds["I_target"].sum(dim="ray")   # this should be multiplied by the area

    return Q_list

def mag(s):
    """return the magnitude of s"""
    return np.sqrt(np.dot(s,s))

def norm(s):
    """return the unit vector in the direction of *s*"""
    return s/mag(s)


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
    ehat = np.zeros((n_dim,n_dim))
    ehat[0] = norm(e_1)
    ehat[1] = norm(e_2)

    if method == "cross":
        ehat[2] = np.cross(ehat[0],ehat[1])
        ehat[1] = np.cross(ehat[2],ehat[0])
    else:
        ehat[1] = ehat[1] - ehat[0]
        ehat[2] = np.cross(ehat[0],ehat[1])

    for i in range(n_dim):
        ehat[i] = norm(ehat[i])

    return ehat


def test():

    x_exit = 207.910 * 1e-3
    x_min = 0.25

    jet_axis = np.array((1,0,0))


    x_center = x_min

    offset = np.array([-0.05,0.0,-0.1])
    pos_center = np.array([x_center,0,0])
    pos_source = pos_center + offset
    pos_target = pos_center - offset

    print("center",pos_center)
    print("source",pos_source)
    print("offset",offset)

    L_beam = mag(offset)*2.0

    beam_axis = -offset
    shat = coordinate_tensor(beam_axis, jet_axis)

    #data_path = r"C:\Users\Huckaby\Desktop\MHD\Simulations\viz_example"
    #fname = "ds_calib.cdf"
    #fname = os.path.join(data_path,fname)
    #f_axi = read_source(fname)

    viz_path = r"C:\Users\Huckaby\Desktop\MHD\Simulations\viz_example"
    print(os.listdir(viz_path))
    fname = "frontCyl_plasma.vtk"
    fname = "frontCyl.vtk"
    fname = os.path.join(viz_path, fname)


    x_beam = np.linspace(0,L_beam)
    d_beam = 1e-3

    tri_out = tri_mesh.new_cylinder(x=x_beam,r=d_beam/2.0,n_theta=60)
    #d = tri_mesh.cylinder(x=x_beam,r=d_beam/2.0)

    ds = tri_out["ds"]
    mesh = tri_out["mesh"]
    #f_axi = create_interpolator(fname, do_plot=False)
    f_axi = new_soln_interpolator(fname, do_plot=False)

    ds["pos_beam"] = ds["pos"]*1.0
    ds["pos"][:] = np.matmul(ds["pos_beam"].to_numpy(),shat) + pos_source
    f_out = f_axi(ds["pos"])
    ds["T"] = ("s","ray"), f_out[:,:,0].T
    ds["K"] = ("s","ray"), f_out[:,:,1].T
    adapt_beam(ds, ["T","K"])

    fig, ax = plt.subplots(2,1,sharex=True)
    ax[0].plot(ds["pos"].isel(ray=5)[:,0],ds["T"].isel(ray=5))
    ax[1].plot(ds["pos"].isel(ray=10)[:,0],ds["K"].isel(ray=10))


    """
    mesh.point_data["points_beam"] = mesh.points*1.0
    new_points = np.matmul(mesh.points,shat) + pos_source
    mesh.points[:] = new_points
    f_out = f_axi(mesh.points)
    mesh.point_data["T"] = f_out[:,0]
    mesh.point_data["K"] = f_out[:,1]
    mesh.point_data["kappa"]
    """

    mesh.save("beam.vtk")



