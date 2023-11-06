"""

Summary
-------

Test several different approaches to generate beam points/volumes
in the plane perpendicular to the beam direction.

1. Rectangular extrusion => hex cells
2. Circular extrusion of in cylindrical coordinate => warped hex cells
3. Triangulate circle and extrude => wedge cells

Different methods are stored in a multiblock vtu file.

The triangulation method, 2 additional beams are created, slightly perturbing
the axial and vertical position of source

Datasets

    cart_beam         cartesian mesh - beam coordinate using 3 position arrays
    cartArray_beam     - beam coordinates using arrays of position vectors
    cart_<method>      - global coordinates using transformation <method>
    cyl_beam          cylindrical mesh - beam coordinates
    cyl_global         - global coodinates

    tri_beam          extruded triangular mesh - beam coordinates
    tri_points         - points without mesh

    tri_x<X>_y<Y>     triangular meshes with different center location
                        <X> - [mm] axial    distance from jet exit
                        <Y> - [mm] vertical offset from jet centerlin

TODO:
    make into function and/or class

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import triangle

from ray_integration import AxiInterpolator
n_dim = 3

def mag(s):
    """return the magnitude of s"""
    return np.sqrt(np.dot(s,s))

def norm(s):
    """return the unit vector in the direction of *s*"""
    return s/mag(s)

def split_dims(x):
    """convert 3D grid of points (4D array) in to a tuple of 3D arrays

    may be good generalize to arbibtrary dimensions
    u = x.T
    v = [u[i].T for i in n_dim]
    which would preserve the ordering

    """
    return x[:,:,:,0], x[:,:,:,1], x[:,:,:,2]

def struct_mesh(x):
    """create structured grid from the 3D grid point points (4D array)

    Notes
    -----
    StructuredGrid class requires that the x,y,z arrays are sent
    """

def cart_transform(x0, shat, x_in, method="matmul"):
    """transform *x_in* to "global" coordinate system


    x_out[i] <=  x0[i] + \sum_j x_in[j] shat[j,i]

    Parameters
    ----------
    x0 : float array [n_dim]
        origin
    shat : float array [n_dim][n_dim]
        shat[i] - unit vector in i-th "local" direction in
         terms of the "global" coordinate vectors
    x_in : 4D float array [*n_s[n_dim],n_dim]
        3D grid of positions

    method : str
        algorithm to use to calculate transformation

    Returns
    -------
    x_out : 4D float array
        3D grid of positions in new coordinate system

    """

    if method=="sum":
        x_out = np.sum(shat[np.newaxis,np.newaxis,np.newaxis,:,:]
                *x_in[:,:,:,:,np.newaxis],axis=-2)
    elif method == "loop_matmul":
        x_out = x_in*0.0
        for i in range(n_dim):
            x_out[:,:,:,i] = np.matmul(x_in[:,:,:,:],shat[:,i])
    elif method == "matmul":
        x_out = np.matmul(x_in[:,:,:,:],shat[:,:])
    elif method == "loop":
        x_out = x_in*0.0
        for i in range(n_dim):
            for j in range(n_dim):
                x_out[:,:,:,j] += xs[:,:,:,i]*shat[i,j]


    return x0 + x_out

def flat_transform(x0, shat, x_in, method="matmul", verbose=False):
    """transfer for array of position vectors, *x_in* using the
    coordinate system defined by, *x0* and *shat*
    """
    if (verbose): print(x_source.shape, shat.shape, x_in.shape)
    x_out = np.matmul(x_in[:,:],shat[:,:])
    if (verbose): print("out:",x_out.shape)
    return x0 + x_out




#  of beam
n_source_pos = 4
n_sub_beams = 3
n_beam = n_source_pos * n_sub_beams

#
# beam perturbation to get  thelocal jet width and spreading rate
#
sub_dy = 3.0e-3
sub_dx = 5.0e-3

# approximate exit
x_exit = 0.21
x_min = 0.25
x_center = np.zeros((n_beam, n_dim))
x_center[0::3,0] = x_min + np.linspace(0.0,0.3,n_source_pos)
x_center[1::3,0] = x_center[0::3,0]
x_center[2::3,0] = x_center[0::3,0] + sub_dx
x_center[1::3,1] = sub_dy

# position of the source relative to the jet cross-point
offset = np.array([-0.05,0.0,-0.1])
x_source = x_center + offset
x_target = x_center - offset

L_beam = mag(offset)*2.0
d_beam = 0.2*L_beam
d_beam = 2e-3
r_beam = d_beam/2.0

jet_axis = np.array((1,0,0))

#
# create a coordinate system centered at the beam source
# oriented to the target
#
shat = np.zeros((n_dim,n_dim))
# beam axis
shat[0] = -offset
#
shat[1] = np.cross(shat[0],jet_axis)
#
shat[2] = np.cross(shat[0],shat[1])

print("Beam unit vectors")
for i in range(n_dim):
    shat[i] = norm(shat[i])
    print("  {:3d} {:5.2f} {:5.2f} {:5.2f}".format(i, *shat[i]))
print()

# number points on the bean axes
n_s = [40,20,35]


#
# cartesian beam
#
s0 = np.linspace(0,1,n_s[0])*L_beam
s1 = np.linspace(-0.5,0.5,n_s[1])*d_beam
s2 = np.linspace(-0.5,0.5,n_s[2])*d_beam

mb = pv.MultiBlock()

# construct using mesh grid
x, y, z = np.meshgrid(s0,s1,s2,indexing='ij')
mesh1 = pv.StructuredGrid(x,y,z)
mb["cart_beam"] = mesh1

xs = np.array([x,y,z])
xs = np.transpose(xs, axes=(1,2,3,0))

# construct using 3D array of position vectors
mesh2 = struct_mesh(xs)
mb["cartArray_beam"] = mesh2

#
# test tranformations
#

transform_methods = "sum loop_matmul matmul loop".split()
i_beam = 0
for method in transform_methods:
    x1 = cart_transform(x_source[i_beam], shat, xs)
    mesh = struct_mesh(x1)
    mb["cart_"+method] = mesh

outpath = "output1"
os.makedirs(outpath, exist_ok="true")
vtk_out = os.path.join(outpath,"transformTests.vtm")
#mb.save(vtk_out)

#
# cylindrical mesh
#
n_r = n_s[1]
n_theta = n_s[2]
r = np.linspace(0,1,n_r)*r_beam
theta = np.linspace(0,1,n_theta)*np.pi*2.0
x_cyl = xs*1.0
x_cyl[:,:,:,1] = r[:,np.newaxis]*np.cos(theta)[np.newaxis,:]
x_cyl[:,:,:,2] = r[:,np.newaxis]*np.sin(theta)[np.newaxis,:]

mb["cyl_beam"] = struct_mesh(x_cyl)
x1 = cart_transform(x_source[i_beam], shat, x_cyl)
mb["cyl_global"] = struct_mesh(x1)

#
# triangular mesh
#

#
# 2D mesh at source
#
n_theta = 30
theta = np.linspace(0,1,n_theta,endpoint=False)*np.pi*2.0
v = np.zeros((n_theta+1,2))
v[1:,0] = np.cos(theta)
v[1:,1] = np.sin(theta)
tri_in = {"vertices": v}
#reg = tri_in["vertices"]*0.0
#reg[:,0] = 1e-3
#reg[:,1] = 1.0
#tri_in["regions"] = reg
n_theta = len(y)
a = np.arange(0,n_theta)
b = a + 1
b[-1] = a[0]
#tri_in["segments"] = np.array([a,b],dtype=int).T
# not needed
plt.plot(v[:,0],v[:,1],'o')

L = 2*np.pi/n_theta
# target area of triangles
area = L**2
# min vertex angle of triangles
angle = 30.0
txt = "a{:f}q{:f}i".format(area, angle)
tri = triangle.triangulate(tri_in, opts=txt)
v = tri["vertices"]
plt.triplot(v[:,0], v[:,1], tri["triangles"])

#
# extrude triangular mesh to wedges
#
n_v, x = v.shape
n_tri, n_verts_per_tri = tri['triangles'].shape

pos = np.zeros((n_s[0],n_v,n_dim))
pos[:,:,0] = s0[:,np.newaxis]
pos[:,:,1] = v[np.newaxis,:,0]*r_beam
pos[:,:,2] = v[np.newaxis,:,1]*r_beam
n_p = n_s[0]*n_v
pos = pos.reshape(n_p,n_dim)
cells = np.zeros((0,0))
celltypes = np.array((0))
#mesh = pv.UnstructuredGrid(cells=cells, celltypes=celltypes, points=pos)

n_x = n_s[0]
cells = np.zeros((n_x-1,n_tri,7),dtype=np.int64)
cells[:,:,1:4] = tri["triangles"]
cells[:,:,4:] = tri["triangles"] + n_v
cells[:,:,0] = 6
cells[:,:,1:] += (np.arange(0,n_x-1)*n_v)[:,np.newaxis,np.newaxis]

n_cells = (n_x-1)*n_tri
n_points = n_x*n_v
celltypes = np.full(n_cells, fill_value=pv.CellType.WEDGE, dtype=np.uint8)
cells = cells.reshape(n_cells,7)
pos = pos.reshape(n_points,n_dim)
mesh = pv.UnstructuredGrid({pv.CellType.WEDGE: cells[:,1:]}, pos)

print("Extruded Triangle Mesh")
print(mesh)
poly = pv.PointSet(pos)
mb["tri_beam"] = mesh
# does not work
mb["tri_points"] = poly

#
# transform to global and create multiple beams
#
results_dir = r"C:\Users\Huckaby\Desktop\MHD\Simulations\viz_example"
sim_run = "."
study = "."
vtk_in = os.path.join(results_dir, study, sim_run, "frontCyl.vtk")
jet = pv.read(vtk_in)

jet_interp = AxiInterpolator(jet)

for i_beam in range(n_beam):
    mesh2 = mesh.copy()
    mesh2.points[:] = flat_transform(x_source[i_beam], shat, mesh2.points)
    data = jet_interp(mesh2.points*1.0)
    n_points, a = mesh2.points.shape


    fmt = "{:5d} {:15s} {:8.1e} {:8.1e}"
    for name, i in jet_interp.index.items():
        if "<" in name:
            #name = '"{}"'.format(name)
            for k in "<>:":
                name = name.replace(k,"_")
            name = name.replace("__","_")

        #print(fmt.format(i, name, data[:,i].min(), data[:,i].max()))
        mesh2.point_data[name] = np.zeros(n_points, dtype=np.float32)
        mesh2.point_data[name][:] = data[:,i]

    x = x_source[i_beam]*1e3
    x[0] -= x_exit*1e3
    mesh_name = "tri_x{:03.0f}_y{:03.0f}".format(x[0],x[1])

    mb[mesh_name] = mesh2
    #mesh2.save(os.path.join(outpath,"{}.vtu".format(mesh_name)))


poly.save(os.path.join(outpath,"triPoints.xyz"))
mb.save(vtk_out)



