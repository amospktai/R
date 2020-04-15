######################################################################
##	Introduction to R for Atmospheric Scientists: Session I
##	Prepared by Amos P. K. Tai (Nov 2011)
##  Updated for Earth System Science students (Apr 2014)
######################################################################

# To execulate ONE command line, Command + Enter (Mac); Ctrl + R (Windows)
# Most important command of all: the HELP function!!!
	?mean

# Get working directory:
	getwd()
	
# Set working directory: (Users: you of course have to set your own working directory; the following is only for my computer...)
	setwd('~/Dropbox/TGABI/R/R_Intro_for_Atmos_Sci/')


######################################################################


# Basic data format and manipulation:

# Vector:
	a = c(21.5, 22.6, 25.4, 27.5, 20.4, 21.8, 24.0)
	a	# Displaying a variable
# Sequence:
	b = seq(1, 7)	# Equivalent to: b = 1:7
	d = seq(0.2, 1.4, by=0.2)
	e = seq(0, 2*pi, length=7)
# Vector of the same value:
	f = rep(1, times=5)
	g = rep(NaN, times=7)	# Empty vector
# Binding vectors to form a matrix:
	A = rbind(a, b, d, e)
	B = cbind(a, b, d, e)
# Matrix:
	C = matrix(1, nrow=3, ncol=5); C
	D = matrix(seq(1, 15), nrow=3, ncol=5); D
	E = matrix(seq(1, 15), nrow=3, ncol=5, byrow=TRUE); E
# Array:
	G = array(0, dim=c(4,3,2))
# Index and assignment:
	D[1,2]
	E[2,]
	G[2,2,2] = 10; G
# Reverse of a vector:
	rev(b)
# Finding which element(s) satisfies a given condition:
	which(b == 1)
	which(b < 3)
	which(b >= 4)
	which(b != 3)
# Listing what you have in the workspace:
	ls()
# Removing variables or functions:
	rm(G)
	rm(list=ls())	# Remove everything!


######################################################################


# Basic mathematical operation:
	1 + 1			# Addition
	2 - 1			# Subtraction
	2*3				# Multiplication
	6/2				# Division
	2^3				# Power

# Note that when a vector or matrix is operated on by another vector or matrix, their dimensions have to match exactly; exception is that a scalar can be operated on any vector or matrix; e.g. b + 1; b + b/2

# Also note that unlike Matlab, multiplication, division and power are by elements. E.g., D*E is not matrix multiplication. Matrix operations are done by:
	H = t(E)		# Tranpose
	J = D %*% H; J	# Matrix multiplication
	D %*% f			# Matrix multiplication
	det(J)			# Determinant
	solve(J)		# Inverse
	eig = eigen(J)	# Eigenvalues and eigenvectors
	eig; eig$values; eig$vectors

# Other basics:	
	sqrt(4)			# Square root
	exp(0.7)		   # Exponential function
	log(2)			# Natural logarithm
	log10(100)		# Base-10 logarithm
	sin(30/180*pi)	# Sine
	cos(60/180*pi)	# Cosine
	tan(45/180*pi)	# Tangent
	asin(0.5)		# Arc-sine
	acos(0.5)		# Arc-cosine
	atan(1)			# Arc-tangent
	factorial(4)	# Factorial
	7.292e-5		# Scientific notation
	6.371e3
	round(3.14)		# Round
	ceiling(3.14)	# Round up
	floor(3.14)		# Round down
	signif(3.14, 2)	# Significant figures
	

######################################################################
	
	
# Basic statistics:

# Loading data:
	load(file='PM25_T_daily_US_201007.RData')
	dim(T.US); head(T.US)
	lon.US; lat.US; day

# Average over all data points:
	mean(T.US, na.rm=TRUE)
# Average over 3rd dimension (day) only to get time average for each lon/lat:
	apply(T.US, c(1,2), mean, na.rm=TRUE)
# Average over 1st and 2nd dimensions (lon/lat) to get regional average for each day:
	apply(T.US, 3, mean, na.rm=TRUE)

# Other important statistics:
	sd(T.US, na.rm=TRUE)		# Standard deviation
	var(T.US, na.rm=TRUE)		# Variance
	median(T.US, na.rm=TRUE)	# Median
	sum(T.US, na.rm=TRUE)		# Sum
	max(T.US, na.rm=TRUE)		# Maximum
	min(T.US, na.rm=TRUE)		# Minimum
	range(T.US, na.rm=TRUE)		# Range
	quantile(T.US, c(0.25, 0.75), na.rm=TRUE)	# Quantile range

# Covariance and correlation:
# E.g. between Boston and San Francisco temperatures:
	cov(T.US[23,9,], T.US[2,7,], use='complete.obs')
	cor(T.US[23,9,], T.US[2,7,], use='complete.obs')
# E.g. between Boston and New York temperatures:
	x = T.US[23,9,]		# Time series in Boston
	y = T.US[21,8,]		# Time series in New York
	cov(x, y, use='complete.obs')
	cor(x, y, use='complete.obs')


######################################################################


# Basic plotting:
	plot(x, y, type='p', xlab='Boston temperature (K)', ylab='New York temperature (K)')
	
# Multiple plots in one single panel:
	T_range = range(c(x, y), na.rm=TRUE)
	plot(day, x, type='l', col='blue', xlab='Date', ylab='Temperature (K)', ylim=T_range)
	matplot(day, y, type='l', col='red', add=TRUE)	# Add another plot on top

# Size of plot region to begin with:
	quartz(title='Temperature', width=4, height=3)		# For Mac OS
	windows(title='Temperature', width=4, height=3)		# For Windows
# You can have multiple panels by setting 'mfrow', e.g., mfrow=c(2,3) will give you 6 panels (2 rows and 3 columns)

# You can customize your plot, e.g. limits, sizes, etc. easily by adding additional parameter, e.g. 'xlim', 'ylim', 'pch', 'lwd', etc. Check out what you can do here:
?par


######################################################################
	
	
# Importing and writing data in ascii or text formats:
	PM25.met.data = read.table(file='PM25_met_Ohio_201007.txt', header=TRUE, sep='\t')
	head(PM25.met.data)		# Showing the first few rows
	tail(PM25.met.data)		# Showing the last few rows
	write.table(PM25.met.data, file='PM25_met_Ohio_201007_new.txt', quote=FALSE, sep='\t')


######################################################################


# Regression analysis:

	PM25 = PM25.met.data$PM25
	Ta = PM25.met.data$Ta
	prec = PM25.met.data$prec
	SLP = PM25.met.data$SLP
	
# Simple linear regression:
	slr = lm(PM25 ~ Ta)		# Regressing PM25 = beta0 + beta1*Ta
	summary(slr)	# Get summary statistics
	plot(Ta, PM25, type='p', xlab='Air temperature (K)', ylab=expression(PM[2.5]*' concentration ('*mu*g*' '*m^-3*')'))
	abline(slr)		# Add regression line; can also add vertical and horizontal lines by setting 'v' and 'h'

# Multiple linear regression:
	mlr = lm(PM25 ~ Ta + prec + SLP)	# Regress PM25 = beta0 + beta1*Ta + beta2*prec + beta3*SLP
	# To exclude intercept: lm(PM25 ~ -1 + T + prec + SLP)
	summary(mlr)

# Extracting statistical info from regression analysis:
	summary(mlr)$coefficients
	summary(mlr)$adj.r.squared
	anova(mlr)

# NOTE: There are many other statistical models in R, e.g. "step" for stepwise regression, "princomp" for principal component analysis, "glm" for generalized linear models, "lmodel2" for reduced major-axis regression, etc. Many of these are not available in built-in packages but can be easily downloaded and installed.


######################################################################


# Saving data or objects into RData format that can be loaded in future:
	save(list=c('PM25', 'Ta', 'prec', 'SLP', 'mlr'), file='PM25_met_regress_result.RData')


######################################################################
##  May leave the following to Session II if time doesn't allow
######################################################################


# There are many packages with useful functions that can be downloaded, installed and loaded. Most useful packages to install: "lmodel2", "ncdf", "fields", "maps", "abind", "mnormt". Suppose you have already installed them, to load them you do:
	library(fields); library(maps)
# E.g. plotting a map for US temperature on 2010-07-01:
	image.plot(lon.US, lat.US, T.US[,,1])
	map('world', add=TRUE)	

# *** For GEOS-Chem users, the most amazing thing is that in IDL there's a routine "bpch2nc" that can convert GEOS-Chem outputs into ncdf format, and R has a package "ncdf" that very efficiently handles ncdf files. I basically do all my GEOS-Chem data analysis in R this way and never really used IDL.


######################################################################


# Topics in Session II:
# Installing and using packages
# Basic programming
# Writing your own functions
# Monte Carlo simulations
# Using lists
# Topics in GIS


######################################################################
######################################################################
