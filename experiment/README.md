# Experiment analysis folder

1. Notebook: analysis of individual experiment dates
currently contains data processing used to create datasets that are collected together ( in `data`) for the final analysis. The data processing needs to be futher refactored and probably just handled along with the data consolidation, or in the main PostProcessing repository. There is some processing that is already happning from the Post Processing repository (run after experiments) for the absorption emission data and all time series data. 
    * In the case of absorption emission (absem) additional processing is appliedm mainly to 2023-05-24 and 2023-05-18 to add the multiplexer (beam location) coordinate and handle that these data were taken with multiple acquisitions between led switching evernts. for all dates this processing is in `absem_setup.py`. This processing happens after the munging in the Post Processing repository and uploading to sharepoint, which mainly reduces the size the data by reducing the wavelength spacing of off-peak wavelength data.  
    * For the microwave scattering data, the `munge_raw_data.py` script processes the raw oscilliscope data (on the NAS Data drive). This mainly reduces the size by reducing the time spacing of points, and aligns data to the correct clock time. 

The analysis scripts are from the initial analysis of the experimental data, and were used for the development of functions (refactored from original code) to collect and add coordinates to data that are used in the final analysis. So the analysis coodes contain extra plots but are now generally obsolete. 

The time signal data (e.g. mass flows, temperatures, etc.) are read from the sharepoint processed data (output of Post Processing repository)

2. data

Collect processed data from the experiment date (notebook) folders. 

Coordinates are added based on the time-series data of relevant signals (mass flow, equivalence ratio, K mass fraction, laser power, and motorized stage position). 

3. analysis 

analayis of final data set collected in `data`

The scripts can be run in an ipython console in VS Code, but can also be rendered to jupyter notebooks with the `render_all.sh` shell script. The notebooks are not included current in the git history. 