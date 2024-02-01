# Uncertainty notes

## definitions 

* (measureable) quantity
    * e.g. length
* value
    *  1.01 cm
* true value -  perfect measurement
* measurement - set of operation to determine the value of a quantity
* ...
* measurand - particular quantity subject to measurement
    * e.g. vapor pressure
* influence quantity - not measurand but affects the result
    * e.g. Temp in lenght measurmeent, freq. in potential measurement
* result - value attributed to measurand, obtained by measurement
* uncorrected result - before correction of systematic error
* accuracy- closeness between result and true value
    * not 'precision'
* repeatability - closeness between results of sucessive measurements with same conditions
* reproducibility - closeness between results of successive measurements with different conditions
    * e.g. different observer, measurement method, ... 


* standard deviation: 


# 2023-12-05

Trying to figure out how to calculate errors after rounds of averaging over mnum and run

* weighted average
* standard deviation or standard deviation of mean?
* lmfit outputs 'standard error' which I believe is sdom

https://blogs.sas.com/content/iml/2019/10/09/statistic-error-bars-mean.html

https://www.statology.org/standard-deviation-vs-standard-error/

```
When to Use Standard Deviation vs. Standard Error

If we are simply interested in measuring how spread out values are in a dataset, we can use the standard deviation.

However, if we’re interested in quantifying the uncertainty around an estimate of the mean, we can use the standard error of the mean.

```


https://lmfit.github.io/lmfit-py/fitting.html

```
After a fit using the leastsq() or least_squares() method has completed successfully, standard errors for the fitted variables and correlations between pairs of fitted variables are automatically calculated from the covariance matrix.

```


questions: 

* are sequential weighted means accurate? I don't think the math works out. 
* is the error formula for weighted average  


