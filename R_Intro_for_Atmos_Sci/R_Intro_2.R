######################################################################
##	Introduction to R for Atmospheric Scientists: Session II
##	Prepared by Amos P. K. Tai (Nov 2011)
##  Updated for Earth System Science students (Apr 2014)
######################################################################


# There are many packages with useful functions that can be downloaded, installed and loaded. Most useful packages to install: "lmodel2", "ncdf4", "fields", "maps", "spam", "abind", "mnormt".

# To install them, either do it via the program menu ("Packages & Data" for Mac version), or use command line:
	install.packages(c('spam', 'fields', 'maps', 'abind', 'ncdf4'))

# Suppose you have already installed them, to load them you do:
	library(fields); library(maps); library(abind)
	
# "fields" and "maps" help you plot maps for gridded data. E.g., we can plot a map for US temperature on 2010-07-01 (see below).

# 'abind' allows you to combine arrays of multiple dimensions along whichever dimension you want.

# *** For GEOS-Chem users, the most amazing thing is that in IDL there's a routine "bpch2nc" that can convert GEOS-Chem outputs into ncdf format, and R has a package "ncdf4" that very efficiently handles ncdf files. I basically do all my GEOS-Chem data analysis in R this way and never really used IDL.

# *** Another very useful package is "ncdf4", which has functions that can read and write NetCDF data format, a widely used format in meteorological and climatological data. E.g., you can download assimilated meteorological data in nc format worldwide from: http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html

# *** An example for how to use "ncdf4" to read, extract and create nc files can be found in "ncdf_hdf5_example.R", also included in this introductory package.

######################################################################


# Plotting maps:
	library(fields); library(maps)
	
# Set working directory:
	setwd('~/Dropbox/TGABI/R/R_Intro_for_Atmos_Sci/')
	
# Loading data:
	load(file='PM25_T_daily_US_201007.RData')
	dim(T.US); head(T.US)
	lon.US; lat.US; day

# Plotting a map for US temperature on 2010-07-01:
	ind = which(day == 20100701)
	image.plot(lon.US, lat.US, T.US[,,ind])
	map('world', add=TRUE)
	
# I've written some efficient functions to plot beautiful maps in the R file "get_geo.R". E.g., "plot.field" and "rwb.colors".
	source('get_geo.R')
	plot.field(T.US[,,ind], lon.US, lat.US, type='def', zlim=c(270, 310))
	

######################################################################
######################################################################


# Basic programming:

# Basically it's about extensively using "for" loops and "if/elseif/else" condition statements


# Example 1:

# Basics using the Ohio PM2.5 example:
	PM25.met.data = read.table(file='PM25_met_Ohio_201007.txt', header=TRUE, sep='\t')
	day = PM25.met.data$day
	PM25 = PM25.met.data$PM25
	Ta = PM25.met.data$Ta
	prec = PM25.met.data$prec
	SLP = PM25.met.data$SLP
	
# Creating an index vector that indicates whether PM2.5 is high, medium or low:
	PM_index = rep(NaN, times=length(PM25))
	for (t in 1:length(PM25)) {
		if (PM25[t] >= 30) PM_index[t] = 'High' else if (PM25[t] >= 10 & PM25[t] < 30) PM_index[t] = 'Med' else PM_index[t] = 'Low'
	}
	print(cbind(day, PM25, PM_index), quote=FALSE)


# Example 2:

# The best way is to illustrate through a more complicated example:
	load(file='PM25_T_daily_US_201007.RData')
	
# Now you want to calculate the regression coefficient for PM25 on T for every 2x2.5 grid square over the US, discarding it if the it is statistically insignificant (p > 0.05) or if there aren't enough data (n < 10):
	reg.coef.PM25.T = matrix(NaN, nrow=length(lon.US), ncol=length(lat.US))
	for (i in 1:length(lon.US)) {
		for (j in 1:length(lat.US)) {
			data.mat = cbind(T.US[i,j,], PM25.US[i,j,])
			data.mat = na.omit(data.mat)	# Manually remove missing data
			if (nrow(data.mat) < 10) {
				reg.coef.PM25.T[i,j] = NaN
			} else {
				X = data.mat[,1]	# T.US
				Y = data.mat[,2]	# PM25.US
				reg = lm(Y ~ X)
				pval = summary(reg)$coefficients[2,4]
				if (pval < 0.05) {
					reg.coef.PM25.T[i,j] = summary(reg)$coefficients[2,1]
				} else {
					reg.coef.PM25.T[i,j] = NaN
				}
			}
		}
	}
	
# You can evaluate multiple condition statements in the following structure:
# if (statement 1) {
	# Do your stuff 1
# } else if (statement 2) {
	# Do your stuff 2
# } else if (statement 3) {
	# etc. etc.
# } else {
	# Do your stuff when all conditions above fail
# }
	
# Now you also want to plot the resulting regression coefficients:
	library(fields); library(maps)
	image.plot(lon.US, lat.US, reg.coef.PM25.T)
	map('world', add=TRUE)
# Or using my written function "plot.field":
	plot.field(reg.coef.PM25.T, lon.US, lat.US, type='sign')

# I have written more customizable and friendly functions to plot beautiful maps in file "get_geo.R", e.g. function "plot.field". You may explore it on your own.

# Saving your precious regression results:
	save(list=c('reg.coef.PM25.T', 'lon.US', 'lat.US'), file='reg_coef_PM25_T.RData')


######################################################################


# Writing your own functions in R:

# Example: finding moving averages for a vector of time series data.
	mov.avg.1D = function(x, t.minus, t.plus) {
		x.ma = rep(NaN, times=length(x))
		for (t in 1:length(x)) {
			t1 = t - t.minus
			t2 = t + t.plus
			if (t1 < 1) t1 = 1
			if (t2 > length(x)) t2 = length(x)
			x.ma[t] = mean(x[t1:t2], na.rm=TRUE)
		}
		return(x.ma)
	}

# Using the written function and plotting the result:
	Ta.ma = mov.avg.1D(Ta, 2, 2)	# 5-day moving averages
	plot(day, Ta, type='o', col='black', xlab='Day of July 2010', ylab='Temperature (K)')
	matplot(day, Ta.ma, type='l', col='red', add=TRUE)

	
# NOTE: I have written many functions to efficiently analyze and plot meteorological, air quality and other multi-dimensional datasets, and also mathematical analysis. They are organized in files: "get_geo.R", "get_met.R", "get_stat.R". To use these functions, simply source the files:
	source(file='get_geo.R')
	source(file='get_stat.R')
	source(file='get_met.R')

# The function "source" basically reads in all commands in a .R file, so you can simply write up an entire script and then source it. 
	
			
######################################################################


# Monte Carlo simulations:

# Basically it refers to any computational algorithm that rely on repeated random sampling to obtain numerical results, typically in order to obtain the uncertainty range (i.e. probability distribution) of an unknown quantity that cannot be otherwise analytically evaluated. R is very efficient in doing such kind of thing.

# For instance, say, A = f(x, y, z), so uncertainty of A should be some combination of that of x, y and x. If A = a*x + b*y + c*z, and if x, y, and z are all normally distributed, A is also normally distributed with variance V(A) = a^2*V(x) + b^2*V(y) + c^2*V(z). However, not all functions are linear, and not all distributions are normal... What if A = exp(x*y - log(z))/(sin(x + y) - tan(z)) + x*y*z, and each of x, y, and z have a different form of distribution? If we have sampling populations for x, y, and z, or if their distributions are known, then we can indeed compuate the distribution of A by repeated sampling of x, y and z.

# But we'll illustrate it by using a linear example just to prove to you that Monta Carlo works, even though we do know the answer analytically.

# Example: Predicting PM25 in Ohio from temperature and precipitation forecast

# From regression, we know that: PM25 = beta0 + beta1*Ta + beta2*prec
# beta0 = -515.165 +/- 15.652 (mean +/- std dev)
# beta1 = 1.7985 +/- 0.0527
# beta2 = -0.4797 +/- 0.0218
# Let's assume beta's are normally distributed (though strictly speaking they're student-t distributed)

# Given tomorrow's weather forecast is such that Ta = 303 K (30 degC), and prec = 2.0 mm/day, what is the uncertainty for PM25 prediction?

# Sampling 10000 times:
	n.trials = 10000
	Ta.predict = 303
	prec.predict = 2.0
	PM25.predict = rep(NaN, times=n.trials)
	for (n in 1:n.trials) {
		# Sampling beta0 assuming normal distribution:
		beta0 = rnorm(1, mean=-515.165, sd=15.652)
		# Sampling beta1:
		beta1 = rnorm(1, mean=1.7985, sd=0.0527)
		# Sampling beta2:
		beta2 = rnorm(1, mean=-0.4797, sd=0.0218)
		# Predicting PM25:
		PM25.predict[n] = beta0 + beta1*Ta.predict + beta2*prec.predict
	}

# Results:
	mean(PM25.predict); sd(PM25.predict)
	hist(PM25.predict)

# Actual mean (known analytically):
	-515.165 + 1.7985*Ta.predict - 0.4797*prec.predict
# Actual standard dev (known analytically):
	sqrt(15.652^2 + 0.0527^2*Ta.predict^2 + 0.0218^2*prec.predict^2)

# If the distribution is unknown but you have a sampling population, you can use function "sample" to sample from the population.


######################################################################


# Root-finding and optimization:

# One-dimensional: e.g. function to find roots or optimize: y = 3*x^2 + 2*x - 5 (we know roots are x1 = 1.0 and x2 = -1.667, and minimum occurs at x = -0.333)
	fx = function(x) 3*x^2 + 2*x - 5
	uniroot(fx, interval=c(0, 2))	# positive root only
	optimize(fx, interval=c(-2, 2))

# We can do multidimensional optimization using function "optim" or through writing our own algorithm:
# E.g., z = x^2 + 2*y^2 - 3*x + 4*y - 5 (we know minimum occurs at x = 1.5, y = -1.0)
	f.xy = function(X) {	# X is a vector of independent variables
		x = X[1]
		y = X[2]
		z = x^2 + 2*y^2 - 3*x + 4*y - 5
		return(z)
	}
	optim(par=c(0, 0), fn=f.xy)
		

######################################################################


# Handling character strings:
	A = 'Hello'
	B = 'World'
	C = paste(A, B, '!', sep=''); C
	C = paste0(A, B, '!'); C
	C = paste(A, B, '!', sep=' '); C
	C = paste0(A, ' ', B, '!'); C
	substr(A, 1, 4)
	strtrim(A, 4)
	paste(substr(A, 1, 3), substr(B, 5, 5), sep='')
# Very useful when you have many different data files to read and thus many filenames to handle.


# Interacting with the operating system (e.g. UNIX in Mac):
	system('ls')

# E.g., duplicating all existing get_xxx.R files and renamed them to be find_xxx.R:
	list.of.files = dir(pattern=c('get', '.R'))
	for (n in 1:length(list.of.files)) {
		command = paste('cp ', list.of.files[n], ' ', 'find', substr(list.of.files[n], 4, 100), sep='')
		system(command)
	}


######################################################################
######################################################################
