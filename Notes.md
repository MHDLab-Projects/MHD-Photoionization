# Notes

## 2023-10-20 

TODO Ideas

* Get 2D data sets integrated
* comparison of cfd and experimental position dependence
    * understand pyvista/paraview
* reorganize/simplify cantera modeling 


Beam measurements 

~ 5mm. Need more accurate measurement


Measuring the mobile absorption system posts with respect to the microwave horns. (distance on torch axis, distance perpindicular). At least 1/8 inch uncertainty. 

Collector (left side): (1 7/8", 4 1/4")
Transmitter (right):  (1 1/2", 3 3/4")

1 + 7/8 = 1.875 
1.875/4.25 = 0.4412 

angle = atan(4.25/1.875)


# 2023-10-27 meeting

Look at K density of all species throughout the length of the jet

Threshold L for experimental nK. Is there an 'effective' L we can determine for different locations to calculate nK. 

look at new interp function from Dave

Make fork of pyvista xarray

# 2023-10-30 meeting

Caught error in AIAA 2020 manuscript. Equation 1 should be 

$ t = \frac{I_{L+F}- I_F}{I_L} = 1 - \alpha $

Linearity Assumption
$ I_{L+F} = I_L + I_F$ 


TODO: just try profiles along the beam centerline to start


# 2023-11-03

Thinking through how to handle experiment notebook

Current setup: process data for each date in experiment notebook, then this repository pulls together all processed data

The processing of data should be ultimately included in the publication repository. 

Option 1: move dates from experiment notebook to publication repository
Option 2: Make a new multi-date processing script. The notebook scripts are heavily refactored in this case and functions used in multi-date

# 2023-11-08 

Chose option 1 