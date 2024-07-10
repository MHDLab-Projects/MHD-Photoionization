* cfd - Analysis of cfd simulations
* viability - PI viability analysis
* em-sim - Mie scattering analyis. 
    - FDTD simulations: meep package: https://meep.readthedocs.io/en/latest/
    - Analytical Mie scattering expressions: PyMieScatt: https://pymiescatt.readthedocs.io/en/latest/
    - For now this is a submodule reference to a separate repository and just generating the figures separately in a conda environment in WSL (see submodule readme) and copying the `meep_sphere/output` folder which contains FDTD datasets and figures used in the SI. This is because the FDTD simulations require conda environment. TODO: reexamine this and automate or document. The FDTD simulations are just a double check for the correct use of analytical mie scattering expressions. They could be removed to automate figure generation in one repository. 