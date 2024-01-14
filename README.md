# Test for proportional harzard ratio hypothesis for python library `statsmodels`

The function is the equivalent of R function survival::cox.zph

You can attached this function to the class `statsmodels.duration.hazard_regression.PHRegResults` and use it as a method.

Normally, you would hope the hypothesis is not rejected so that you can assume the effect of the variable(s) in the Cox regression is not time-variant.

