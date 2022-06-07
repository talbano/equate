# equate_2.0.8 (2022-06-05)

* Added sample method for freqtab object, useful for bootstrapping.  

* Converted vignette to R markdown.  

* Small updates throughout.  

# equate_2.0.7 (2018-04-06)

* Update to summary method for bootstrap output, to correct the calculation of mean standard error, bias, and root mean square error over score points.  

# equate_2.0-6 (2017-01-09)

* bootstrap function can now extract presmoothing arguments, if present, to run presmoothing once for all functions requesting it. Note that, when requested, presmoothed score distributions will be used for all methods.  

* New fit statistics for loglinear presmoothing model comparison: likelihood ratio chi-square, Freeman-Tukey chi-square, consistent AIC, Cressie-Read, and Goodman.  

* Fix plot axes in plot.freqtab where scales do not match observed score range.  

* Legend text fixed for single-group and equivalent-groups designs in plot.equate and plot.bootstrap.  

* Algorithm for default model comparison in loglinear presmoothing was fixed to address problem with degrees of 10 or larger.  

* Bug in Levine true score equating, introduced in 2.0.5, was fixed.  

# equate_2.0-5 (2016-10-20)

* Support for linear equating in a single-group design with counterbalancing.  

* Analytic standard errors added for linear equating under single-group counterbalanced and equivalent-groups designs, and chained, Tucker, and Levine methods with nonequivalent-groups.  

* Analytic standard errors added for equipercentile equating under equivalent-groups design.  

* Fix bug in equipercentile equating for scores with zero frequencies.  

* Missing bracket fixed in loglinear function.  

* Vignette updated to match recent JSS publication.  

# equate_2.0-4 (2016-4-26)

* Support for equating scales that have inconsistent score increments.  

* summary methods for equate and bootstrap objects fixed to support partial output, e.g., only standard errors.  

* The equating design can be specified when creating a frequency table object. See the design argument in the freqtab functions.  

* Another bug in summary method for equate objects fixed. An error was returned for anything but equipercentile functions.  

* Support added for automatic selection of best fitting model in loglinear presmoothing function. See the choose argument in presmoothing and loglinear.  

* Small tweaks and improvements to help files.  

* Development version of package now available at https://github.com/talbano/equate  

* Package is now built using RStudio and roxygen2.  

# equate_2.0-3 (2014-11-1)

* Functionality has been added for equating under a multiple anchor test design. New methods include nominal weights and Tucker multivariate and multivariate frequency estimation for Braun/Holland and equipercentile equating.  

* The freqtab class has been restructured, with a new, more general approach to the construction and manipulation of multivariate frequency tables. A freqtab object is now a simple table array rather than a data.frame.  

* formula method added for loglinear presmoothing via the presmoothing function.  

* Loglinear presmoothing supports multivariate frequency distributions, e.g., distributions with multiple anchor tests. Polynomials and interactions can be specified using a matrix of power terms, in addition to the score function or the degree (now 'degrees') argument.  

* Impossible scores can be removed before loglinear presmoothing using the 'rmimpossible' argument.  

* Model comparisons for loglinear presmoothing are based on the 'models' argument, which is a numeric vector with different levels specifying the sequential models to be fitted. Also, when 'stepup' is TRUE, the default models are now chosen based on common powers across variables.  

* The default legend for plotting bootstrap output was fixed.  

* Summary method for bootstrap output can take a subset of the score scale. The method was also updated to include weighted mean RMSE.  

# equate_2.0-2 (2014-2-21)

* Some minor cleanup.  

* plot.bootstrap was fixed for objects containing output for only one equating.  

* print method was added for bootstrap objects.  

# equate_2.0-1 (2014-2-1)

* points.freqtab fixed to handle plotting frequencies that are all the same.  

* Bug fixed in printing and summarizing equipercentile output.  

* Score conversion with chained equipercentile output now works. Previously, there was an error in how the smoothing method was interpreted.  

* Score conversion with equipercentile now works when x only contains values with percentile ranks of 0 or 1. Previously, this produced an error.  

# equate_2.0-0 (2014-1-7)

* New "general linear" equating type.  

* Single-group equating is now supported, where "x" is bivariate and "y" is missing.  

* Equipercentile equating will accept a rescaled "y" scale, e.g., with non-integer scores. The score increment still must be constant across the scale.  

* The identity linking function is modified using the arguments "lowp" and "highp", which were previously only used in circle-arc equating.  

* "highp" is now used to find "Ky" in equipercentile equating. The "Ky" argument is no longer used.  

* Identity equating is now considered a linear equating type, though it is still specified with type = "identity".  

* New "composite" object class is used to create composite equatings, i.e., weighted combinations of two or more other equating functions. Synthetic equating, combining an identity function with another function, is no longer supported within the equate function.  

* "px" and "fx" are now internal, meaning a help file is no longer included.  

* "yx" is no longer included separately in equating output. Equated scores can be accessed through "concordance$yx" instead.  

* Analytical standard errors are under construction and not currently included in equating output.  

* For loglinear smoothing within the equate function, a single "scorefun" argument has replaced the separate "xscorefun" and "yscorefun" arguments. The equate function cannot be used to smooth "x" and "y" with different score functions.  

* Bug in the calculation of the intercept in chained mean and linear equating was fixed.  

* Circle-arc equating with nonequivalent groups was fixed.  

* Chained equipercentile score conversion was fixed.  

* Score conversion with frequency estimation equipercentile was fixed.  

* "summary" method replaces "descript" for summarizing frequency tables.  

* Additional descriptive statistics are now included in equating output.  

* Levine true score equating printed output was fixed.  

* Bootstrap standard error calculation was fixed for scales of different lengths.  

* "plot.equate" was function modified.  

* "plot.freqtab" was function updated to improve plotting logliner presmoothing results.  

* "freqtab" function was updated to handle missing data.  

# equate_1.2-0 (2013-9-20)
  
* Package vignette updated.  

* freqtab function updated to accept item response data.  

* loglinear function updated - original estimation algorithm has been dropped, and the function is now a wrapper for 'glm'.  

* new function plot.equate, for plotting and comparing equating results.  

* equate and other functions modified to handle linking of forms with different score scales.  

* Bug fixed in the frequency estimation method, where marginal frequencies are zero.  

* Help files updated.  

# equate_1.1-4 (2011-8-23)

* Package vignette updated.  

* Print method fixed.  

* equate function fixed for equipercentile conversion of a vector of scores to the scale defined in a previous equating.  

# equate_1.1-3 (2011-2-10)

* verbose option added to loglinear smoothing.  

* equate function updated to support conversion of a vector of scores to the scale defined in a previous equating.  

# equate_1.1-2 (2011-1-7)

* Circle-arc equating fixed.  

* Package vignette updated.  

# equate_1.1-1 (2010-12-15)

* New equating type: synthetic, a combination of any other equating type and identity equating.  

* Circle-arc equating fixed.  

* skew.freqtab and kurt.freqtab fixed.  

# equate_1.1-0 (2010-11-15)

* plot.freqtab updated.  

* New equating type: nominal weights.  

* loglinear smoothing function fixed.  

# equate_1.0-0 (2010-03-24)

* New package vignette introduces concepts and demonstrates basic functionality.  

* New equating type: circle-arc.  

* New equating methods: Levine true score, Braun/Holland, and chained equating.  

* equate function converts a vector of scores to another scale based on output from a previous equating.  

* Proportional weighting when estimating the synthetic distribution, using argument 'w'.  

* Layout of freqtab function changed.  

* Simple plotting method added for univariate and bivariate frequency tables.  

* se.boot function returns matrix of bootstrap replications.  

* Equipercentile equating updated to handle frequency tables missing substantial amounts of data.  

* Functions added for descriptive statistcs of frequency tables.  

* Rd files updated and consolidated.  

# equate_0.1-1 (2009-11-05)

* Functions loglinear and freqbump extended to handle smoothing of bivariate distributions.  

* equating functions accept scores only as frequency tables (univariate or bivariate).  

* Updated DESCRIPTION and Rd files.  

# equate_0.1-0 (2009-11-01)

* initial release on CRAN.  
