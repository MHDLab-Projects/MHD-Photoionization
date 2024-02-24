"""

Summary
-------
Create surface and extruded volume meshes for cylinders using "triangle",
a 2D triangulation tool


"""
import numpy as np
import matplotlib.pyplot as plt
import triangle
import pyvista as pv
import xarray as xr

def tri_to_pv(tri):
    pass

def tri_to_xr(tri):
    pass

def tri_correct_area(tri):
    """calculate and apply curvature correction to outer triangles"""

    pos = tri["vertices"]
    tri_verts = tri["triangles"]
    radius = np.sum(pos**2,axis=1)
    outer_verts = abs(radius  - 1) < 1e-3
    outer_tris = np.where(np.sum(outer_verts[tri_verts],axis=1) == 2)[0]
    n_outer_tris = outer_tris.shape[0]
    def mag(x):
        return np.sum(x**2, axis=-1)**0.5

    tri["area_corr"] = tri["area"]*0.0

    # list of vertices
    v_out = tri_verts[outer_tris]
    n_outer_tris = len(v_out)
    for i in range(n_outer_tris):
        v_side = v_out[i][ np.where(outer_verts[v_out[i]]) ]
        v_in = v_out[i][ np.where(outer_verts[v_out[i]]== False) ]
        p_side = pos[v_side]
        r = mag(p_side)
        #print(r)
        val = np.dot(p_side[0],p_side[1])/(r[0]*r[1])
        alpha = np.arccos(val)
        A_corr = 0.5*(alpha - np.sin(alpha))*r[0]**2
        #print(i, alpha, val, A_corr, r)

        tri["area_corr"][outer_tris[i]] = A_corr

    tri["area"] += tri["area_corr"]


def tri_calc_area(tri):
    """given triangulation *tri* calculate the area of the triangles

    ref:
    """

    pos = tri["vertices"]
    tri_verts = tri["triangles"]

    p1 = pos[tri_verts[:,0]]
    p2 = pos[tri_verts[:,1]]
    p3 = pos[tri_verts[:,2]]

    pA = (p2 - p1)
    pB = (p3 - p2)
    pC = (p1 - p3)

    sA = np.sum(pA**2,axis=-1)**0.5
    sB = np.sum(pB**2,axis=-1)**0.5
    sC = np.sum(pC**2,axis=-1)**0.5
    s = (sA + sB + sC)/2.0
    A = (s*(s-sA)*(s-sB)*(s-sC))**0.5
    tri["area"] = A

def tri_calc_weights(tri, verbose=False):
    """given triangulation *tri* calculate the area of the cell around
    a vertex and weights"""

    if not "area" in tri:
        tri_calc_area(tri)

    n_v = len(tri["vertices"])
    w = np.zeros(n_v)
    n_cells = np.zeros(n_v,dtype="int")
    n_tris, n_vert_per_tri = tri["triangles"].shape
    #print("n_verts per triange = ",n_vert_per_tri)
    #print("# node cells", n_cells.min(), n_cells.max() )

    A = tri["area"]
    for i, vert_i in enumerate(tri["triangles"]):
        w[vert_i] += A[i]/3.0

    for i, vert_i in enumerate(tri["triangles"]):
        n_cells[vert_i] += 1

    tri["v_area"] = w
    tri["v_weights"] = w/np.sum(w)
    tri["v_triangles"] = n_cells


def tri_plot(tri, ax=None):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    if ax == None:
        fig, ax = plt.subplots(1,3)
        [a.set_aspect(1.0) for a in ax]

    v = tri['vertices']
    cmap = plt.get_cmap("viridis")
    ax[0].triplot(v[:,0], v[:,1], tri["triangles"])
    for i, name in enumerate(["v_area","v_triangles"]):
        w = tri[name]
        divider = make_axes_locatable(ax[i+1])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        out = ax[i+1].scatter(v[:,0], v[:,1], c=w, cmap=cmap)
        plt.colorbar(out, cax=cax)
    fig.tight_layout()
    return fig, ax


def new_poly_disk(n_r=3, n_sides=6, outer_points="double", correct_area=True, do_plot=False):
    """create triangular mesh in/on unit circle/disk based on regular polygon point spacing

    Parameters
    ----------
    n_r : int
        number of radial cells
    n_sides : int
        number of sides of the polygon
    outer_points : string
        spacing algorithm for the point on the circle
        "double" - for the ring use twice the number points as the next
            inner ring

    correct_area : bool
        apply curvature correction with calculating the area of the
        outer triangles

    Returns
    -------
    tri : dictionary
        output of triangle plus
        area[n_tri]       : triangle area
        vert_area[n_v]    :
        vert_weights[n_v] :

    """
    n_dim2 = 2
    pos = np.zeros((1,n_dim2))
    pos_list = [pos]
    for i_r in range(1,n_r+1):
        n_theta = n_sides*i_r
        r = float(i_r)/n_r
        n_theta = n_sides*i_r
        theta = np.linspace(0,1.0,n_theta,endpoint=False)*np.pi*2.0
        if (outer_points == "double") and (i_r == n_r):
            n_theta = 2*n_sides*(i_r-1)
            theta = np.linspace(0,1.0,n_theta,endpoint=False)*np.pi*2.0
            #theta += (theta[1]-theta[0])*0.25

        #print(i_r, n_theta)

        pos = np.zeros((n_theta,n_dim2))
        pos[:,0] = r*np.cos(theta)
        pos[:,1] = r*np.sin(theta)
        pos_list.append(pos)

    pos_all = np.concatenate(pos_list)
    tri_in = {"vertices": pos_all}
    tri = triangle.triangulate(tri_in)
    tri_calc_area(tri)

    # correct area
    if correct_area:
        tri_correct_area(tri)

    tri_calc_weights(tri)

    if do_plot:
        tri_plot(tri)

    return tri



def new_disk(n_theta=30, angle=30.0, area_scale=1.1, correct_area=True, do_plot=False):
    """create triangle mesh on unit disk

    Parameters
    ----------
    n_theta : int, optional
        DESCRIPTION. The default is 30.
    angle : float, optional
        DESCRIPTION. The default is 30.0.
    area_scale : float
        scaling factor for triangle target area

    Returns
    -------
    None.

    """
    theta = np.linspace(0,1,n_theta,endpoint=False)*np.pi*2.0
    v = np.zeros((n_theta+1,2))
    v[1:,0] = np.cos(theta)
    v[1:,1] = np.sin(theta)
    tri_in = {"vertices": v}

    L = 2*np.pi/n_theta
    # target area of triangles
    area = area_scale*3.0**0.5/2.0*L**2

    tri_opts = "a{:f}q{:f}iv".format(area, angle)
    tri = triangle.triangulate(tri_in, opts=tri_opts)

    tri_calc_area(tri)
    if correct_area:
        tri_correct_area(tri)
    tri_calc_weights(tri)


    #print("# node cells", n_cells.min(), n_cells.max() )
    #print("# node areas", A.min(), A.max() )
    #print("# node weights", w.min(), w.max() )

    if do_plot:
        tri_plot(tri)

    return tri


n_dim = 3
def new_cylinder(x=[0,1], r=1.0, n_theta=30, angle=30.0, area_scale=1.1, correct_area=True, do_plot=False):
    tri = new_disk(n_theta, angle, area_scale, correct_area, do_plot)
    tri['vertices'] *= r
    tri['area'] *= r**2
    tri['v_area'] *= r**2
    d = new_extrusion(tri, x, do_plot)
    d.update({"tri":tri})
    return d


def tri_print_center_area(tri):
    t = tri['triangles']
    v0 = 0
    center_tris = np.where(np.any(t==v0, axis=1))
    print("\ncenter data")
    print("vertex", tri["vertices"][v0] )
    print("tris", tri["triangles"][center_tris] )
    print("area",tri["area"][center_tris] )
    print()

def new_poly_cylinder(x=[0,1], r=1.0, n_r=3, n_sides=6, outer_points="double", correct_area=True, do_plot=False):
    tri = new_poly_disk(n_r=n_r, n_sides=n_sides, outer_points=outer_points,
                        correct_area=correct_area, do_plot=do_plot)
    tri['vertices'] *= r
    tri['area'] *= r**2
    tri['v_area'] *= r**2
    d = new_extrusion(tri, x, do_plot)
    d.update({"tri":tri})
    return d


def new_extrusion(tri, x=[0,1], do_plot=True):
    v = tri['vertices']
    n_x = len(x)
    n_v, n_d = v.shape
    pos = np.zeros((n_x, n_v, n_dim))

    pos[:,:,0] = x[:,np.newaxis]
    pos[:,:,1] = v[np.newaxis,:,0]
    pos[:,:,2] = v[np.newaxis,:,1]

    # transform to unstructured mesh
    n_tri, n_verts_per_tri = tri['triangles'].shape
    cells = np.zeros((n_x-1,n_tri,7),dtype=np.int64)
    cells[:,:,1:4] = tri["triangles"]
    cells[:,:,4:] = tri["triangles"] + n_v
    cells[:,:,0] = 6
    cells[:,:,1:] += (np.arange(0,n_x-1)*n_v)[:,np.newaxis,np.newaxis]

    n_cells = (n_x-1)*n_tri
    n_points = n_x*n_v
    #celltypes = np.full(n_cells, fill_value=pv.CellType.WEDGE, dtype=np.uint8)
    cells = cells.reshape(n_cells,7)
    pos1 = pos.copy().reshape(n_points,n_dim)
    mesh = pv.UnstructuredGrid({pv.CellType.WEDGE: cells[:,1:]}, pos1)

    #print(pos1.min(axis=1),pos1.max(axis=1))

    coords = {"s":x,"ray":np.arange(n_v),"i_dir":"x y z".split()}
    ds = xr.Dataset(coords)
    ds["pos"] = ("s","ray","i_dir"), pos*1.0

    ds["v_area"] = "ray", tri["v_area"]
    ds["v_weights"] = "ray", tri["v_weights"]
    ds["v_triangles"] = "ray", tri["v_triangles"]

    # this is view not a copy
    for i, name in enumerate("xyz"):
        ds.coords["{}".format(name)] = ("s","ray"), ds["pos"].data[:,:,i]

    n_points, n_tmp = mesh.points.shape
    u = np.zeros((n_x,n_v),'d')
    u[:,:] = ds["v_area"]
    mesh.point_data["v_area"] = u.reshape(n_points)

    u_tri = np.zeros((n_x-1,n_tri),'d')
    u_tri[:,:] = tri["area"][np.newaxis,:]
    mesh.cell_data["tri_area"] = u_tri.reshape(n_cells)

    return {"raw":pos,"ds":ds,"mesh":mesh}






