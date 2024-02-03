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

# 2023-11-27

TODO: add spectral of calibration. 

initial estimate of beam propogation path

error analysis: average led on/off over entire test case. 


# 2023-12-01

Compare LED to expected thorlabs

measure drift over hours. 

calculate error in fit from parameter covariance, and compare with error bars

# 2023-12-07 

## Uncertainty analysis

Sources of uncertiainty 

* variation within a given 'run': intra-run
    * time fluctuation
    * process unsteadiness
    * dimenion 'mnum' (aka measurement number)
* Variation between runs: inter-run
    * Process changes (emulsion etc.)
    * dimension: 'run' 
* fitting
    * fitting covariance

Questions:

How to handle combining errors between inter and intra process?
* appears inter-run variation is larger
* the 'true value' of the distibution is changing
* Weighted average? Only one average/stddev calculation over all measurements?

How to integrate fitting and associated parameter uncertianty?
* Do we perform averages before or after fitting? 3 possibilities. 
* Believe error corresponds to standard error

Is error in weighted average the standard error or standard deviation

$ \sigma_{stderr}  = \sigma_{stdev}/\sqrt{counts}$ 

Believe the counts used in this should be the counts of intra-run



#TODO: look at theoretical curves of alpha, wrt various inlcuding temperature

#TODO test effect of data point spacing on absem fit: does data point number correspond to counts?

Add a known noise 


#TODO: compare fitting grouping
1. average intra-run before fit
2. fit each individual acquistion then average fit parameters
3. 'global' least squares fitting. fit all acquisition to the same curve with one standard error 

# 2024-01-02

## Pipeline packages

Looking around for a simple python package for defining data processing/analysis pipelines. Seems like most are designed for specific types of data and more complex than needed. 

* Sklearn
    * https://scikit-learn.org/stable/modules/generated/sklearn.pipeline.Pipeline.html
    * too specific to machine learning 

* Web-based. To much of an emphasis on scheduling...
    * Luigi
        * https://luigi.readthedocs.io/en/stable/
    * prefect
        * https://www.prefect.io/
    * apache airflow
        * https://codesolid.com/airflow-python-etl/

* Bonobo
    * https://www.bonobo-project.org/
    * looks promising, but last commit 2021
    * http://docs.bonobo-project.org/en/master/faq.html#why-not-use-some-library-instead

* found an astrophysics package with internal implementation like what I'm thinking
    * https://pythonhosted.org/Astropysics/coremods/pipeline.html#astropysics.pipeline.Pipeline

* Pipeln
    *https://cgarciae.gitbook.io/pypeln/readme



* build from scractch...
    * https://dkraczkowski.github.io/articles/crafting-data-processing-pipeline/
    * why is this not a package


# 2023-01-08

TODO: 

- [x] Improve handling of nans for global fit. 
- [x] Then develop second mws fitting pipeline and revisit absem pipe 2. 
    * Had to drop all wavelength for absem, which leads to distorted spectrum for low kwt.
    * different length t_eval error message is showing up for mws

- [x]  Need to process multiple values for ne0 for cfd data. 
    * interpolate to kwt values of experiment. 

- [] Add mw_horn nK fit to parameter for comparison. 

- [x] add basic exponential fit

- [] look at OH-


# 2024-01-19

Look at phase oscillations

look at CFD e- Kp and OH

# 2024-02-02

- [x] Main K dependence plot with CFD time windows. Develop functionality for sub timewindows. 

- [] plot of error for each indivisual test case 

- [] Collect physical properties

- [x] check surfactant emulsion recipe 
    use test plan for recipe and add surfactant to CFD input
