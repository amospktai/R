# This file contains useful functions for analyzing time series data, and also meteorological data obtained from the NCEP/NCAR Reanalysis Project.
# Last update: Mar 2017, Amos P. K. Tai (amostai@cuhk.edu.hk)

###############################################################################

# Data extraction from .RData file:
# This function simply extracts an data array "data.array" from a given .RData file associated with "filename".
# Both "data.array" and "filename" are a string scalar. "data.array" is the name given to the variable of interest stored in the .RData file associated with "filename".
load.RData = function(filename, data.array) {
	.env = new.env()
	load(filename, envir=.env)
	data = get(data.array, envir=.env)
	return(data)
}

###############################################################################

# Function extraction from "functions.RData" file:
# This function simply extracts a function "function" from "functions.RData" file in the fixed location depending on "server" which specifies either one of the location: "group" = Jacob Group server; "local" = personal computer.
# "function" is a string scalar, and is the name given to the function of interest stored in the "functions.RData" file.
load.functions = function(func, server="local") {
	.env = new.env()
	if (server == "local") filename = "~/Dropbox/Research/get_data/functions.RData"
	if (server == "group") filename = "~/get_data/functions.RData"
	load(filename, envir=.env)
	f = get(func, envir=.env)
	return(f)
}

###############################################################################

# Monthly statistics (most general form):
# This function calculates the monthly statistics for one year of data, containing exactly 365 or 366 entries.
# Default "fun" is "mean".
monthstats1 <- function(data, fun=mean) {
    if (length(data) == 365) {
    	M = c(fun(data[1:31], na.rm=TRUE), fun(data[32:59], na.rm=TRUE), fun(data[60:90], na.rm=TRUE), fun(data[91:120], na.rm=TRUE), fun(data[121:151], na.rm=TRUE), fun(data[152:181], na.rm=TRUE), fun(data[182:212], na.rm=TRUE), fun(data[213:243], na.rm=TRUE), fun(data[244:273], na.rm=TRUE), fun(data[274:304], na.rm=TRUE), fun(data[305:334], na.rm=TRUE), fun(data[335:365], na.rm=TRUE))
    } else {
        M = c(fun(data[1:31], na.rm=TRUE), fun(data[32:60], na.rm=TRUE), fun(data[61:91], na.rm=TRUE), fun(data[92:121], na.rm=TRUE), fun(data[122:152], na.rm=TRUE), fun(data[153:182], na.rm=TRUE), fun(data[183:213], na.rm=TRUE), fun(data[214:244], na.rm=TRUE), fun(data[245:274], na.rm=TRUE), fun(data[275:305], na.rm=TRUE), fun(data[306:335], na.rm=TRUE), fun(data[336:366], na.rm=TRUE))
    }
    return(M)
}

###############################################################################

# Monthly statistics for multidimensional data in x/y/time format:
# This function calculates the monthly statistics for one year of data for multiple locations, and utilizes another function "monthstats1". X is an array with dimensions of: dim1 = x-coord (e.g. longitude); dim2 = y-coord (e.g. latitude); dim3 = time (in days).
# Default "fun" is "mean".
monthstats2 = function(X, fun=mean) {
    monthly = array(0, dim=c(dim(X)[1:2], 12))
    for (i in 1:(dim(X)[1])) {
        for (j in 1:(dim(X)[2])) {
            monthly[i,j,] = monthstats1(data=X[i,j,], fun=fun)
        }
    }
    return(monthly)
}

###############################################################################

# Annual statistics from daily data:
# This function calculates the annual statistics from start year "year1" to end year "year2" from a vector of daily data "X". There must be no missing days in X from year1 to year2.
yearstats1 = function(X, year1, year2, fun) {
	year = seq(year1, year2, by=1)
	annual = rep(0, times=length(year))
	if ((year1/4 - floor(year1/4)) == 0) {
		day.pat = rep(c(366, 365, 365, 365), times=1000)
	} else if ((year1/4 - floor(year1/4)) == 0.75) {
		day.pat = rep(c(365, 366, 365, 365), times=1000)
	} else if ((year1/4 - floor(year1/4)) == 0.5) {
		day.pat = rep(c(365, 365, 366, 365), times=1000)
	} else {
		day.pat = rep(c(365, 365, 365, 366), times=1000)
	}
	annual[1] = fun(X[1:day.pat[1]], na.rm=TRUE)
	day.cum = cumsum(day.pat)
	for (t in 2:length(annual)) annual[t] = fun(X[(day.cum[t-1]+1):day.cum[t]], na.rm=TRUE)
	return(annual)
}

###############################################################################

# Calculate central moving averages for time series data:
# This function calculates the central moving averages for each point of a vector of time series data.
# "time.low" and "time.up" both have absolute values and have the same units with "time", and thus (time.low + time.up) defines the moving averaging window.
mov.avg = function(data, time, time.low, time.up) {
	m.avg = rep(0, times=length(data))
	for (t in 1:length(time)) {
		ind = which(time >= (time[t]-time.low) & time <= (time[t]+time.up))
		m.avg[t] = mean(data[ind], na.rm=TRUE)
		}
	return(m.avg)
}

###############################################################################

# Calculate moving averages for 3-D time series data:
# This function calculates the moving averages for each point of 3-D time series data with time being the 3rd dimension.
# "dk1" and "dk2" define the window for averaging: e.g. dk1 = 15 and dk2 = 15 find the moving average over 15 time points before to 15 time points after a given time point.
# This function does not allow non-uniform time intervals.
mov.avg.3D = function(data, dk1, dk2) {
	m.avg = array(NaN, dim=dim(data))
	for (t in 1:dim(data)[3]) {
		k1 = t - dk1
		k2 = t + dk2
		if (k1 < 1) k1 = 1
		if (k2 > dim(data)[3]) k2 = dim(data)[3]
		m.avg[,,t] = apply(data[,,k1:k2], c(1,2), mean, na.rm=TRUE)
	}
	return(m.avg)
}

###############################################################################

# Calculate weighted, central moving averages for time series data:
# This function calculates the weighted, central moving averages for each point of a vector of time series data.
# "weights" is a vector defining the weights given to the data in the averaging window, e.g. c(1, 2, 3, 2, 1) gives the center point of the averaging window a weight of 3, the data points right before and after the center a weight of 2, and the data points +2 and -2 points away from the center a weight of 1, etc.
# The length of "weights" necessarily also defines the averaging window.
# The time series cannot contain NaN values and has to have uniform intervals.
mov.avg.wgt = function(data, weights) {
	dum = length(weights)/2 - 0.5
	data.long = c(rep(NaN, times=dum), data, rep(NaN, times=dum))
	m.avg = rep(0, times=length(data))
	for (t in 1:length(data)) {
		sum.weights = sum(weights*(data.long[t:(t+dum*2)]/data.long[t:(t+dum*2)]), na.rm=TRUE)
		m.avg[t] = sum(data.long[t:(t+dum*2)]*weights/sum.weights, na.rm=TRUE)
	}
	return(m.avg)
}

###############################################################################

# Calculate wind direction from wind velocity vectors:
# This function calculates the wind direction in angle ("theta", in radians) starting counterclockwise from the east, from wind velocity components U and V.
find.wind.dir = function(U, V) {
	theta = atan(V/U)
	L1 = (U > 0 & V > 0)
	L2 = (U < 0 & V > 0)
	L3 = (U < 0 & V < 0)
	L4 = (U > 0 & V < 0)
	L5 = (U > 0 & V == 0)
	L6 = (U == 0 & V > 0)
	L7 = (U < 0 & V == 0)
	L8 = (U == 0 & V < 0)
	if (length(which((L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8) > 1)) > 0) stop('Sum of Ls greater than 1!')
	theta = theta*L1 + (theta + pi)*L2 + (theta + pi)*L3 + (theta + 2*pi)*L4 + 0*pi*L5 + 0.5*pi*L6 + 1*pi*L7 + 1.5*pi*L8
	# Obselete code:
    #if (U > 0 & V > 0) theta <- theta
    #else if (U < 0 & V > 0) theta <- theta + pi
    #else if (U < 0 & V < 0) theta <- theta + pi
    #else if (U > 0 & V < 0) theta <- theta + 2 * pi
    #else {}
	return(theta)
}

###############################################################################

# Convert between NCEP/NCAP (Julian) time format (hour since year 0000) and more common YYYYMMDD format:
# These two functions convert between the two formats, and uses a reference year 1948 with 1948/01/01 = 17067072.
# Please note that the date is strictly for 00:00-00:00 UTC. Make the necessary time zone shift for converting to and from local date.

from.yyyymmdd = function(date) {
	# "date" and "time" cannot be a vector or matrix, and cannot be in the year of 1948.
	year = round(date/10000)
	month = round((date-year*10000)/100)
	day = date-year*10000-month*100
	year.diff = year-1948
	year.ind = seq(1, year.diff)
	add.day = 0
	for (i in 1:length(year.ind)) {
		if (((i-1)/4) == round((i-1)/4)) add.day = add.day + 366
		else add.day = add.day + 365
	}
	month.ind = seq(1,month)
	for (j in 1:length(month.ind)) {
		if (j == 1) add.day = add.day + 0
		else if (j == 2 | j == 4 | j == 6 | j == 8 | j == 9 | j == 11) add.day = add.day + 31
		else if (j == 5 | j == 7 | j == 10 | j == 12) add.day = add.day + 30
		else if (j == 3 & (year.diff/4) == round(year.diff/4)) add.day = add.day + 29
		else add.day = add.day + 28
	}
	add.day = add.day + day
	time = 17067072 + (add.day-1)*24
	return(time) 
}

to.yyyymmdd = function(time) {
	# "date" and "time" cannot be a vector or matrix, and cannot be in the year of 1948.
	add.day = floor((time-17067072)/24) + 1
	year4 = ceiling(add.day/(366+365*3)) - 1
	year = 1948
	if (year4 == 0)	{
		add.day = add.day
		year = year
	} else {
		add.day = add.day - year4*(366+365*3)
		year = year + year4*4
	}
	y = c(366, 365, 365, 365)
	y = cumsum(y)
	if (add.day <= y[1]) {
		add.day = add.day
		year = year
	}
	else if (add.day > y[1] & add.day <= y[2]) {
		add.day = add.day - y[1]
		year = year + 1
	}
	else if (add.day > y[2] & add.day <= y[3]) {
		add.day = add.day - y[2]
		year = year + 2
	}
	else {
		add.day = add.day - y[3]
		year = year + 3
	}
	if (((year-1948)/4) == round((year-1948)/4)) m = c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
	else m = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
	m = cumsum(m)
	if (add.day <= m[1]) {
		add.day = add.day
		month = 1
	}
	else if (add.day > m[1] & add.day <= m[2]) {
		add.day = add.day - m[1]
		month = 2
	}
	else if (add.day > m[2] & add.day <= m[3]) {
		add.day = add.day - m[2]
		month = 3
	}
	else if (add.day > m[3] & add.day <= m[4]) {
		add.day = add.day - m[3]
		month = 4
	}
	else if (add.day > m[4] & add.day <= m[5]) {
		add.day = add.day - m[4]
		month = 5
	}
	else if (add.day > m[5] & add.day <= m[6]) {
		add.day = add.day - m[5]
		month = 6
	}
	else if (add.day > m[6] & add.day <= m[7]) {
		add.day = add.day - m[6]
		month = 7
	}
	else if (add.day > m[7] & add.day <= m[8]) {
		add.day = add.day - m[7]
		month = 8
	}
	else if (add.day > m[8] & add.day <= m[9]) {
		add.day = add.day - m[8]
		month = 9
	}
	else if (add.day > m[9] & add.day <= m[10]) {
		add.day = add.day - m[9]
		month = 10
	}
	else if (add.day > m[10] & add.day <= m[11]) {
		add.day = add.day - m[10]
		month = 11
	}
	else {
		add.day = add.day - m[11]
		month = 12
	}
	day = add.day
	date = year*10000 + month*100 + day
	return(date)
}

###############################################################################

# Make a vector of dates:
make.date.vec = function(start.date, end.date, return="date.vec") {
	# This function creates a vector of dates specified by "start.date" and "end.date", both of which are numeric, e.g. 20040101 for Jan 1, 2004.
	# This function utilizes other functions "from.yyyymmdd" and "to.yyyymmdd".
	time.vec = seq(from.yyyymmdd(start.date), from.yyyymmdd(end.date), by=24)
	date.vec = time.vec
	for (t in 1:length(date.vec)) date.vec[t] = to.yyyymmdd(time.vec[t])
	if (return == "date.vec") return(date.vec)
	if (return == "time.vec") return(time.vec)
}

###############################################################################

# Match two datasets according to dates:
match.date = function(Y, Y.date, X, X.date) {
	# This function finds shared dates between two time-series datasets "Y" and "X", eliminates unshared dates, and returns a subset of original data sharing common dates.
	# "Y" is a numeric vector; "X" can be a numeric vector or a matrix with time/date being the 1st dimension (by row).
	# "Y.date" and "X.date" are the dates corresponding to "Y" and "X" respectively.
	# Output includes: "$date" are the shared dates; "$Y" and "$X" are the subsets of original data sharing common dates.
	date.vec = sort(unique(c(X.date, Y.date)))
	if (is.vector(X)) MAT = matrix(NaN, nrow=length(date.vec), ncol=2)
	else MAT = matrix(NaN, nrow=length(date.vec), ncol=(1+ncol(X)))
	for (t in 1:length(date.vec)) {
		Y.ind = which(Y.date == date.vec[t])
		if (length(Y.ind) == 0) MAT[t,1] = NaN
		else MAT[t,1] = Y[Y.ind]
		X.ind = which(X.date == date.vec[t])
		if (length(X.ind) == 0) {
			if (is.vector(X)) MAT[t,2] = NaN
			else MAT[t,2:(1+ncol(X))] = NaN
		}
		else {
			if (is.vector(X)) MAT[t,2] = X[X.ind]
			else MAT[t,2:(1+ncol(X))] = X[X.ind,]
		}
	}
	MAT = cbind(date.vec, MAT)
	MAT = na.omit(MAT)
	date = MAT[,1]
	Y = MAT[,2]
	X = MAT[,-c(1,2)]
	OUTPUT = list(date=date, Y=Y, X=X)
	return(OUTPUT)
}

###############################################################################

# Day of year from date and vice versa:
# This function finds the day of year given a standard date yyyymmdd, or finds the date in yyyymmdd from the day of year given a year.
# Input argument "yyyymmdd", "day.of.year" or "yyyy" can be a vector.

date.to.day = function(yyyymmdd, leap=FALSE) {
	if (leap) days.in.month = c(0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) else days.in.month = c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
	day = yyyymmdd - signif(yyyymmdd, 6)
	month = (yyyymmdd - signif(yyyymmdd, 4) - day)/100
	day.of.year = day + cumsum(days.in.month)[month]
	return(day.of.year)
}

day.to.date = function(day.of.year, yyyy, leap=FALSE) {
	if (leap) days.in.month = c(0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) else days.in.month = c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
	month = NULL
	for (n in 1:length(day.of.year)) month = c(month, max(which(cumsum(days.in.month) < day.of.year[n])))
	day = day.of.year - cumsum(days.in.month)[month]
	yyyymmdd = yyyy*1e4 + month*1e2 + day
	return(yyyymmdd)
}

###############################################################################

# Leap year:
# This function determine whether a given year is a leap year:
is.leap = function(yyyy) {
	if (yyyy/4 == round(yyyy/4)) leap = TRUE else leap = FALSE
	return(leap)
}
	
###############################################################################

# Seasonal data extraction:
# This function extracts seasonal data from "spdata", which is either a 1-D time-series numeric vector, 2-D data matrix where time/date is the 1st dimension (row), or a 3-D data array where time/date is the 3rd dimension.
# "date.vec" is a numeric vector specifying the dates corresponding to the 1-D time series, 1st (for 2-D) or 3rd (for 3-D) dimension of "spdata" in the form of YYYYMMDD (numeric).
# "season" is a string specifying the season of interest, either one of "MAM" (spring), "JJA" (summer), "SON" (fall) or "DJF" (winder).

find.season = function(spdata, date.vec, season) {
	date.vec.new = date.vec - signif(date.vec, digits=4)
	if (season == "MAM") ind.season = which(date.vec.new >= 301 & date.vec.new <= 531) else if (season == "JJA") ind.season = which(date.vec.new >= 601 & date.vec.new <= 831) else if (season == "SON") ind.season = which(date.vec.new >= 901 & date.vec.new <= 1130) else if (season == "DJF") ind.season = which(date.vec.new >= 1201 | date.vec.new <= 229) else ind.season = 1:length(date.vec)
	if (is.vector(spdata)) spdata.new = spdata[ind.season] else if (length(dim(spdata)) == 2) spdata.new = spdata[ind.season,] else spdata.new = spdata[,,ind.season]
	return(spdata.new)
}

###############################################################################

# Seasonal mean statistics:
# This function finds the seasonal mean values for "spdata", a 1-D time-series numeric vector, 2-D data matrix where time/date is the 1st dimension (row), or a 3-D data array where time/date is the 3rd dimension.
# "date.vec" is a numeric vector specifying the dates corresponding to the 1-D time-series, 1st (for 2-D) or 3rd (for 3-D) dimension of "spdata" in the form of YYYYMMDD (numeric).
# "season" is a string specifying the season of interest, either one of "MAM" (spring), "JJA" (summer), "SON" (fall) or "DJF" (winder).

find.season.mean = function(spdata, date.vec, season) {
	date.vec.new = date.vec - signif(date.vec, digits=4)
	if (season == "MAM") ind.season = which(date.vec.new >= 301 & date.vec.new <= 531) else if (season == "JJA") ind.season = which(date.vec.new >= 601 & date.vec.new <= 831) else if (season == "SON") ind.season = which(date.vec.new >= 901 & date.vec.new <= 1130) else if (season == "DJF") ind.season = which(date.vec.new >= 1201 | date.vec.new <= 229) else ind.season = 1:length(date.vec)
	if (is.vector(spdata)) season.mean = mean(spdata[ind.season], na.rm=TRUE) else if (length(dim(spdata)) == 2) season.mean = apply(spdata[ind.season,], 2, mean, na.rm=TRUE) else season.mean = apply(spdata[,,ind.season], c(1,2), mean, na.rm=TRUE)
	return(season.mean)
}

###############################################################################

# Monthly data extraction:
# This function extracts monthly data from a 3-D data array "spdata", where time or date is the 3rd dimension.
# "date.vec" is a numeric vector specifying the dates corresponding for the 3rd dimension of "spdata" in the form of YYYYMMDD (numeric).
# "month" is a numeric specifying the month of interest, e.g. for Jan month = 1 and for Dec month = 12.

find.month = function(spdata, date.vec, month) {
	date.vec.new = date.vec - signif(date.vec, digits=4)
	ind.month = which(date.vec.new >= (month*100 + 1) & date.vec.new <= (month*100 + 31))
	spdata.new = spdata[,,ind.month]
	return(spdata.new)
}

###############################################################################

# Monthly mean statistics:
# This function finds the monthly mean values for a 3-D data array "spdata", where time or date is the 3rd dimension.
# "date.vec" is a numeric vector specifying the dates corresponding for the 3rd dimension of "spdata" in the form of YYYYMMDD (numeric).
# "month" is a numeric specifying the month of interest, e.g. for Jan month = 1 and for Dec month = 12.

find.month.mean = function(spdata, date.vec, month) {
	date.vec.new = date.vec - signif(date.vec, digits=4)
	ind.month = which(date.vec.new >= (month*100 + 1) & date.vec.new <= (month*100 + 31))
	monthly.mean = apply(spdata[,,ind.month], c(1,2), mean, na.rm=T)
	return(monthly.mean)
}

###############################################################################

# Saturation water vapor pressure:
# This function calculates the saturation water vapor pressure in [hPa] for a given T or a vector/matrix of T in [K] based on Lowe and Ficke (1974).
satvap.H2O = function(T, give.derivative=FALSE) {
	T = T - 273.15			# Convert to deg C.
	a0 = 6.107799961
	a1 = 4.436518521e-1
	a2 = 1.428945805e-2
	a3 = 2.650648471e-4
	a4 = 3.031240396e-6
	a5 = 2.034080948e-8
	a6 = 6.136820929e-11
	es = a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4 + a5*T^5 + a6*T^6
	if (give.derivative) {
		des.dT = a1 + 2*a2*T + 3*a3*T^2 + 4*a4*T^3 + 5*a5*T^4 + 6*a6*T^5
		return(des.dT)
	} else {
		return(es)
	}
}

###############################################################################

# Coriolis parameter:
# This function calculates the Coriolis parameter in [s^-1] for a given latitude in [deg].
f_Cor = function(lat) 2*7.292115e-5*sin(lat/180*pi)

###############################################################################

# Plot two variables with different units in the same plot:
plot.x2 = function(X, Y1, Y2, type1='l', type2='l', col1='blue', col2='red', xlab='', ylab1='', ylab2='', ylab2.line=1.4, axis2=NULL, axis2.ticks=5) {
	Y1.scaled = (Y1 - sum(range(Y1, na.rm=TRUE))/2)/diff(range(Y1, na.rm=TRUE))
	Y2.scaled = (Y2 - sum(range(Y2, na.rm=TRUE))/2)/diff(range(Y2, na.rm=TRUE))
	plot(X, Y1, type=type1, col=col1, xlab=xlab, ylab=ylab1)
	Y2.plot = Y2.scaled*diff(range(Y1, na.rm=TRUE)) + sum(range(Y1, na.rm=TRUE))/2
	matplot(X, Y2.plot, type=type2, col=col2, add=TRUE)
	if (is.null(axis2)) {
		Y2.axis = signif(seq(signif(range(Y2, na.rm=TRUE)[1], digits=3), signif(range(Y2, na.rm=TRUE)[2], digits=3), length=axis2.ticks), digits=3)
	} else {
		Y2.axis = axis2
	}
	Y2.axis.scaled = (Y2.axis - sum(range(Y2, na.rm=TRUE))/2)/diff(range(Y2, na.rm=TRUE))
	Y2.axis.plot = Y2.axis.scaled*diff(range(Y1, na.rm=TRUE)) + sum(range(Y1, na.rm=TRUE))/2
	axis(side=4, at=Y2.axis.plot, labels=Y2.axis)
	mtext(text=ylab2, side=4, line=ylab2.line)
}

###############################################################################

# Find time-lagged correlations between two variables:
# "x" and "y" are two vectors of variables for which the time-lagged correlations are to be found.
# "n" is the maximum number of time lags considered, e.g. if n = 5, correlations will be found for x(t) vs. y(t), x(t) vs. y(t-1), x(t) vs. y(t-2), ..., x(t) vs. y(t-5).
cor.lag = function(x, y, n, plot.cor=FALSE) {
	cor.vec = rep(NaN, times=n+1)
	for (N in 1:(n+1)) {
		x.subset = x[N:length(x)]
		y.subset = y[1:(length(y)-N+1)]
		cor.vec[N] = cor(x.subset, y.subset, use='complete.obs')
	}
	if (plot.cor) plot(0:n, cor.vec, 'o', xlab='time lag', ylab='correlation coefficient')
	return(cor.vec)
}

###############################################################################

# Convert between GCM or CTM hybrid eta levels, pressure levels and altitudes:
# This function converts between hybrid eta levels, pressure levels and altitudes that are commonly used to GCM or CTM, e.g. GEOS, GEOS-Chem.
# For now, only "GEOS-5" vertical grid is represented.
# "lev.in" is the input vector of levels, with type specified by "from"; "to" is the type you want "var.in" to be converted to.
# "from" and "to" can be "eta", "pres" and "alt".
level.convert = function(lev.in, from='eta', to='pres', model='GEOS-5') {
	if (model == 'GEOS-5') {
		eta.frame = c(0.000000, 0.000028, 0.000055, 0.000127, 0.000199, 0.000399, 0.000599, 0.001109, 0.001619, 0.002816, 0.004013, 0.006588, 0.009162, 0.014342, 0.019523, 0.023755, 0.027987, 0.033814, 0.039641, 0.047641, 0.055641, 0.066559, 0.077477, 0.084313, 0.091149, 0.099191, 0.107233, 0.116695, 0.126157, 0.137287, 0.148418, 0.161513, 0.174608, 0.190061, 0.205513, 0.223772, 0.242032, 0.263587, 0.285142, 0.309854, 0.334566, 0.353349, 0.372133, 0.390927, 0.409720, 0.428528, 0.447337, 0.466153, 0.484970, 0.503795, 0.522620, 0.541449, 0.560278, 0.579115, 0.597953, 0.616790, 0.635628, 0.654471, 0.673314, 0.685878, 0.698442, 0.711006, 0.723570, 0.736134, 0.748698, 0.761265, 0.773832, 0.786400, 0.798967, 0.809021, 0.819075, 0.826616, 0.834157, 0.841698, 0.849239, 0.856781, 0.864323, 0.871864, 0.879406, 0.886948, 0.894489, 0.902031, 0.909573, 0.917116, 0.924658, 0.932200, 0.939743, 0.947285, 0.954828, 0.962370, 0.969913, 0.977456, 0.984999, 0.992500, 1.000000)
		pres.frame = c(0.010, 0.038, 0.066, 0.139, 0.211, 0.414, 0.617, 1.134, 1.651, 2.864, 4.077, 6.685, 9.293, 14.542, 19.792, 24.080, 28.368, 34.272, 40.175, 48.282, 56.388, 67.450, 78.512, 85.439, 92.366, 100.514, 108.663, 118.250, 127.837, 139.115, 150.393, 163.661, 176.930, 192.587, 208.244, 226.745, 245.246, 267.087, 288.927, 313.966, 339.005, 358.038, 377.070, 396.112, 415.155, 434.212, 453.269, 472.335, 491.401, 510.475, 529.550, 548.628, 567.706, 586.793, 605.880, 624.967, 644.054, 663.146, 682.239, 694.969, 707.699, 720.429, 733.160, 745.890, 758.621, 771.354, 784.088, 796.822, 809.556, 819.743, 829.929, 837.570, 845.211, 852.852, 860.493, 868.135, 875.776, 883.418, 891.059, 898.701, 906.342, 913.984, 921.626, 929.268, 936.911, 944.553, 952.195, 959.837, 967.480, 975.122, 982.765, 990.408, 998.051, 1005.650, 1013.250)
		alt.frame = c(80.581, 72.180, 68.392, 63.053, 59.924, 54.834, 51.788, 47.135, 44.286, 40.166, 37.574, 34.024, 31.716, 28.654, 26.596, 25.307, 24.240, 23.020, 22.004, 20.836, 19.855, 18.727, 17.773, 17.243, 16.753, 16.222, 15.731, 15.198, 14.706, 14.170, 13.674, 13.134, 12.633, 12.086, 11.578, 11.021, 10.504, 9.936, 9.409, 8.846, 8.320, 7.943, 7.582, 7.237, 6.905, 6.585, 6.277, 5.980, 5.692, 5.413, 5.142, 4.879, 4.623, 4.375, 4.132, 3.896, 3.665, 3.439, 3.219, 3.074, 2.932, 2.792, 2.654, 2.517, 2.382, 2.249, 2.118, 1.988, 1.860, 1.759, 1.659, 1.584, 1.510, 1.436, 1.363, 1.290, 1.218, 1.146, 1.075, 1.004, 0.934, 0.864, 0.795, 0.726, 0.657, 0.589, 0.521, 0.454, 0.387, 0.320, 0.254, 0.189, 0.123, 0.058, -0.006)
		if (from == 'eta') x = eta.frame
		if (from == 'pres') x = pres.frame
		if (from == 'alt') x = alt.frame
		if (to == 'eta') y = eta.frame
		if (to == 'pres') y = pres.frame
		if (to == 'alt') y = alt.frame
		lev.out = spline(x, y, xout=lev.in)$y
	} else {
		stop('For now other model grids are not supported.')
	}
	return(lev.out)
}
	
###############################################################################

	