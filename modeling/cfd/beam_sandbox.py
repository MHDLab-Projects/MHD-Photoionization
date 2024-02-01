"""
Summary
-------
General testing of pyvista and other tools

"""
import os
import pyvista as pv
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import numpy as np


results_dir = r"C:\Users\Huckaby\Desktop\MHD\Simulations\viz_example"
case = "."
study = "."
fname = os.path.join(results_dir, study,case,"frontCyl.vtk")

mesh = pv.read(fname)
n_dim = 3
x_min = 0.2
x_max = 0.5
dx = 0.05
n_beam = int((x_max-x_min)/dx)

source_offset = np.array([-0.05, 0, -0.1])
source_distance = np.sum(source_offset**2)**0.5
# location where beam goes through center plane (z=0)
beam_cL = np.zeros((n_beam,n_dim))
beam_cL[:,0] = np.linspace(x_min, x_max, n_beam)
beam_cL[:,1] = 0.0

beam_source = beam_cL + source_offset
beam_dir = -source_offset/source_distance
beam_L = 2*source_distance
beam_target = beam_source + beam_dir*beam_L

print("source:",beam_source)
print("target:",beam_target)
print("CL    :",beam_cL)

beam_lines = [mesh.sample_over_line(beam_source[i],beam_target[i],resolution=100) for i in range(n_beam)]

beam_line0 = []
beam_line1 = []
os.makedirs("output1", exist_ok=True)
for i in range(n_beam):
    line = pv.Line(beam_source[i],beam_target[i],2)
    print(line)
    line.save("output1/line_{}.vtk".format(i))
    sp = pv.Spline([beam_source[i],beam_cL[i],beam_target[i]],2)
    print(sp)
    sp.save("output1/spline_{}.vtk".format(i))


    line1 = mesh.slice_along_line(line)
    print(line1)

    beam_line0.append(line)
    beam_line1.append(line1)

    line1.save("output1/beam_{}.vtk".format(i))



