# AgeDLModels
Code for simulating data from and fitting age-specific distributed lag models

This repository contains code to simulate from the age-specific distributed lag models proposed in Heaton et al (2018) Environmetrics with a few notable changes.  First, in order to decrease the number of parameters needs to simulate the data, this code will simulate from a Poisson distribution. Second, we assume that there is no effect for sex or age so no parameters are needed for these effects.  Finally, we do not specify an offset.

To simulate data, source the file SimulateFromAgeDL.R with your choices for the following parameters:
 	- min.age - the youngest aged person in your dataset
	- max.age - the oldest aged person in your dataset
	- L - the maximum non-zero lag
	- M - the maximum zero lag
	- r - the reduction factor specifying the number of basis functions
	- nu - the smoothness parameter for the Gaussian process prior
	- alpha.age - the decay parameter for the Matern correlation function across age
	- alpha.lag - the decay parameter for the matern correlation function across lag
	- sigma2 - the marginal *unconstrained* variance of the DL coefficients
	- beta0 - the intercept for the model
Once these values are specified, the code will simulate (i) DL coefficients from the constrained PP prior then (ii) simulate data y ~ pois(exp(beta0+DL.vals)) where DL.vals is the lagged temperatures times the simulated DL coefficients.  The lagged temperatures are taken as the average daily temperatures from Houston as in the aforementioned paper.

Once data has been simulated, you can fit the simulated data using NIMBLE from FitAgeDLwithNIMBLE.R.  You need only source this file for it to run.