# Workshop: occupancy models in R

A workshop on occupancy models in R. This workshop uses the 'unmarked' package to construct occupancy models in R. All examples use included simulated datasets.

This workshop covers:

* single season occupancy models and goodness-of-fit tests
* building and comparing multiple models
* model-averaging predictions for occupancy models
* fitting dynamic (multi-season) occupancy models

Contents:

* 01_installingoccupancypackages.R checks if required packages are installled. Installs them if they are not.
* 02_singleseasonoccupancy.R contains code for introduction to occupancy models
* 03_singleseasonoccupancy_part2_modelcomparison.R contains code for building and comparing single season occupancy models
* 04_dynamicoccupancy.R contains code for building dynamic (multi-season) occupancy models
* detection_history.csv contains simulated detection history data
* site_cov.csv contains simulated site covariates (forest cover and agriculture cover)
* effort.csv contains simulated observation covariate for search effort per survey
* observers.csv contains simulated observation covariate for number of observers per survey
* dynamic_detection_history.csv contains simulated detection history data used in 04_dynamicoccupancy
* dynamic_site_cov.csv contains simulated site covariates (wetland, temp, prec) used in 04_dynamicoccupancy
* dynamic_effort.csv contains simulated observation covariate for search effort per survey used in 04_dynamicoccupancy
