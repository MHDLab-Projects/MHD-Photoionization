# -*- coding: utf-8 -*-
"""
Summary
-------
pvista utilities for axisymmetric meshes


Classes
-------
AxiMesh - wrapper to perform 3D sampling of an axi-symmetric mesh using
    pyvista sampling tools

pv_interpolator - creates a 3D interpolating function over selected variables
    of a mesh

AxiInterpolator - wrapper to perform 3D sampling of an axi-symmetric mesh using
    scipy interpolation routines

adaptive_line_interpolator - create an interpolator over a line based on a source
    dataset adjusting the number of points based on gradients of the data in
    the source
    
Notes
-----

Links
-----
https://docs.pyvista.org/version/stable/api/core/_autosummary/pyvista.DataSetFilters.sample.html
https://docs.pyvista.org/version/stable/examples/01-filter/resample.html#resampling-example


"""
import os
import numpy as np
import scipy.interpolate as interpolate
from scipy.constants import R as R_SI, N_A as N_A_SI
import pyvista as pv
import matplotlib.pyplot as plt

from mhdpy.pyvista_utils import AxiMesh

# convert from SI (m-mol-s) to MKS (m-kmol-s)
R_u = R_SI * 1e3
N_A = N_A_SI * 1e3

def pv_interpolator(mesh, var_names=[]):
    p = mesh.points
    names = []
    data = []
    if len(var_names) == 0:
        vars = mesh.point_data.keys()
    else:
        vars = [k for k in var_names if (k in mesh.point_data)]

    for k in vars:
        v = mesh.point_data[k]
        if len(v.shape) == 1:
            data.append(v)
            names += [k]
        if len(v.shape) == 2:
            n_points, n = v.shape
            for i in range(n):
                data.append(v[:,i])
                names.append("{}_{}".format(k,i))
    data = np.array(data).T
    return interpolate.LinearNDInterpolator(p[:,:2], data), names
       
        
class AxiInterpolator:

    def __init__(self, mesh, var_names=[]):
        self.mesh = mesh
        self.interp, self.names = pv_interpolator(mesh, var_names)
        self.index = { x : self.names.index(x) for x in self.names }


    def __call__(self, points):

        p = points.T
        x = p[0]
        r = (p[1]**2 + p[2]**2)**0.5

        return self.interp((x,r))
    
    def eval_dict(self, points):        
        u = self.interp((x,r))
        
        return 
    
def line_interpolator(f_pos, f_val, a=0.0, b=1.0, max_iter=10, rtol=0.01, do_plot=False, verbose=1):
    """grid interpolator *f_val* to a line interpolator

    Parameters
    ----------
    f_pos : function(s)
        calculate cartesian coordintes from line coordinate, *s*
    f_val : function(pos[:,n_dim])[:,n_var]
        calculate array of n_var values, 
    a : float
        starting point of line
    b : interval end
        end point of line
    max_iter : int
        number of iterations
    rtol : float
        max relative difference between values
    verbose : bool
    
    do_plot : bool
        plot data
        
    Returns
    -------


    Notes
    -----
    Given:
        1. f_val(x) - a "value" function which calculates a set of values at cartesian
            position, *x*
        2. f_pos(s) - a "position" function which calculates a position given a line
            position, *s*

    1. finds a set of points, s[] on the interval between a and b, such that
        the difference in the values between any two adjacents points on s[]
        is less than a specied value
        
        s[i] - set of n_points in [0,1]
        x[i,i_dim] - cartersiam value
      
     2. adds news points between s[i] and s[i+1] if 
        any( abs(u[i,i_var] - u[i+1,ivar] ) ) > w[ivar]*rtol


    """
    s = np.linspace(a,b,21)
    x = f_pos(s)
    v = f_val(x)
    #
    w_small = 1e-6
    w_vsmall = 1e-12
    
    # weighting function
    w = (v.max(axis=0) - v.min(axis=0)) + (v.max(axis=0) + v.min(axis=0))*w_small + w_vsmall
    
    if do_plot:
        cmap = plt.colormaps["viridis"]
        n_points, n_vars = v.shape
        fig, ax = plt.subplots(n_vars, 2)
        for i in range(n_vars):
            ax[i,0].plot(s,v[:,i],'k.')

    if verbose:
        header = 50*"-"
        print(header)
        print("Adpative Interpolator")
        print("{:4} {:5} {:4} {:>6} {:>9} {:>9}".format("iter","points","#new","%new","maxDelta","maxLoc"))
        print(header)
    for i in range(max_iter):

        dv = abs(v[1:,:] - v[:-1,:])
        ii = np.where( dv > rtol*w[np.newaxis,:] )
        i_new = np.unique(ii[0])
        if verbose:
            n_new = len(i_new)
            n_s = len(s-1)
            f = 100.0*n_new/n_s
            dv_norm = np.max(dv/w,axis=-1)
            print("{:4d} {:5d} {:4d} {:6.1f} {:9.1e} {:9.1e}".format(i, n_s, n_new, f, dv.max(), s[dv_norm.argmax()]))
        if len(ii[0]) < 1:
            break
        s_new = (s[i_new] + s[i_new+1])*0.5
        x_new = f_pos(s_new)
        v_new = f_val(x_new)

        s = np.concatenate([s,s_new])
        v = np.concatenate([v,v_new])

        ii = np.argsort(s)
        v = v[ii]
        s = s[ii]
        dv = abs(v[1:,:] - v[:-1,:])
        if do_plot > 0:
            color = cmap( float(i)/max_iter )
            for i in range(n_vars):
                if do_plot == 1:
                    ax[i,0].plot(s_new[:],v_new[:,i],'.', color=color, alpha=0.05)
                else:
                    ax[i,0].plot(s[:],v[:,i],'-', color=color, alpha=0.05)
                ax[i,1].hist(dv[:,i], color=color, orientation="horizontal")
                ax[i,1].set_ylabel(r"$\Delta v$"+"_{}".format(i))
                ax[i,0].set_ylabel("v_{}".format(i))
                ax[-1,0].set_xlabel("beam coordnate")
                ax[-1,1].set_xlabel("# of new points")



    x = f_pos(s)
    f = interpolate.interp1d(s, v.T)
    if verbose:
        n_s = len(s)
        a = np.arange(n_s)
        for i in range(n_vars):
            ax[i,0].set_xlabel("s")
            tw = ax[i,0].twiny()
            tw.plot(a,v[:,i],'.', alpha=0.02, label="i")
            tw.set_xlabel("i")
        print(header)
        fig.tight_layout()
        return f, fig
    return f

    
if __name__ == "__main__":    
    
    results_dir = r"C:\Users\Huckaby\Desktop\MHD\Simulations\viz_example"
    sim_run = "."
    study = "."
    fname = os.path.join(results_dir, study,sim_run,"frontCyl.vtk")
    n_dim = 3
    mesh = pv.read(fname)

    p = mesh.point_data
   
    p['C_mix'] = p['p']/(p['T']*R_u)
    # TODO calculate mole fraction
    p['C_K'] = p['C_mix']*p['K']
    f_mesh = pv_interpolator(mesh, ["C_K","T"])

    f_axi = AxiInterpolator(mesh, ["C_K", "T"])
    
    
    x_CL = np.zeros((9,n_dim))
    x_CL[:,0] = np.linspace(0.2,0.6,9)
    
    x_offset = np.array([-0.05,0.0,-0.1])
    x_source = x_CL + x_offset
    x_target = x_CL - x_offset
    L = np.sum((x_source - x_target)**2,axis=1)
    
    
    
    
    
    
    