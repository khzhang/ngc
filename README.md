# ngc:
an R Package for estimating network Granger causality. It implements the general lasso, truncating lasso, adaptively thresholded lasso, and group lasso regression regularization framework for longitudinal and time series data. These methods are described in Shojaie and Michailidis 2010, Shojaie and Michailidis 2011 and Shojaie and Michailidis 2015.

# Installation
Run the following lines in R.
```R
install.packages("devtools")
devtools::install_github("khzhang/ngc", build_vignettes=T)
```

# What this package does
Estimates, prediction, and visualizes.
Use `ngc' function to estimate the Granger Causality.
```R
fit <- ngc(X,  method = 'truncate',
            typeIerr = 0.05)
```
Use `predict' function to predict the following time points.

```R
predict(fit, 2)
```
Use `plot.ngc' function to visualize the Granger Causality.
```R
plot.ngc(fit3)
```


# Demo
A vignette is available [here](Vignette/metrics_eval.html). The vignette gives a tutorial about the `ngc' packages and demonstrates the usage of different arguments. In addition, a comparison between the implemented methods is available [here](Vignette/metrics_eval.html). This file shows the comparison between the implemented estimation methods, including average running time and relative errors. To view these files, you can either download files in your current working directory and open in R, or you can open the html from [Chrome](https://www.google.com/chrome/) or [Firefox](https://www.mozilla.org/firefox/). 


# Usage
Install the R package, and in R call the 
```R
ngc::ngc()
```
function on the dataset; an example is given [here](demo/demo.R). The data could be either an array with dimension nXpXT or a matrix with dimension pXT. n is the number of observations, p is the number of concurrent time series, and T is the number of time points to be considered. 

There are three main methods with different options. Users can specify the methods to do the estimation. Check vignette (khzhang/ngc/Vignette/Introduction-to-ngc.html) for more details. 


# References
Shojaie A. and Michailidis G. (2010) Discovering Graphical Granger Causality Using a Truncating Lasso Penalty, Bioinformatics, 26(18): i517-i523

Shojaie A., Basu S. and Michailidis G. (2012) Adaptive Thresholding for Reconstructing Regulatory Networks from Time Course Gene Expression Data, Statistics In Biosciences 4(1): 66-83

Basu S., Shojaie A. and Michailidis G. (2015) Network Granger Causality with Inherent Grouping Structure, Journal of Machine Learning Research (JMLR)

