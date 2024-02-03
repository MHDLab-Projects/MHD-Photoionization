"""
Summary
-------


Routines
--------
beam -
    (In progress) objected oriented Beam class and functions

beam_utils -
    functions to do beam integration using xarray datasets as the data container

jet_profiles -
    analytical jet profile

mesh_gen -
    Test several different approaches to generate beam points/volumes
    in the plane perpendicular to the beam direction.

parameter_estimation -
    Sensitivity analysis and parameter fitting (round-trip) based on an
    analytical jet profile.

pv_axi_utils -
    pvista utilities for axisymmetric meshes

ray_integration (broken) -
    testing scipy tools for beam integration with absorption

tri_mesh -
    Create surface and extruded volume meshes for cylinders using "triangle",
    a 2D triangulation tool

Notebooks
---------

beam_sensitivity
    1. reads cfd data and makes interpolation functions
    2. reads and smooths source wavelength vs. intensity profile
    3. creates beam.Beam and interpolates cfd data to beam w/ adpative refinement of mesh
    4. calculates theoretical absorption cross section
    5. TODO beam calculation

xrBeam_integration
    1. reads cfd data and makes interpolation functions
    2. creates beam.Beam and interpolates cfd data to beam w/ adpative refinement of mesh
    3. reads and smooths source wavelength vs. intensity profile
    4. calculates
        - fictious gaussian absorption cross section
        - voigt absorption cross section
    5. calculates beam intensity and total intensity at target for both cross-section functions


PolyMesh_scratchPad
    - notebook for devoloping algorithms and manipulate meshes

test__tri_mesh
    - testing and demonstration of tri_mesh modules

"""


