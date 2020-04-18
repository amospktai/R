# This file contains useful functions for general numerical and statistical analysis.
# Last update: Dec 2016, Amos P. K. Tai (amostai@cuhk.edu.hk)

################################################################################

# Area under curve:
trapz = function(x, y) {
	# This function calculates the area under a curve y = f(x), i.e. integral of f(x)dx, using linear approximation.
  	A.i = (y[-1]+y[-length(y)])*(x[-1]-x[-length(x)])/2
  	A = sum(A.i)
  	return(A)
}

################################################################################

# Roots of quadratic equation:
quadroot = function(a, b, c) {
	# This function finds the roots of the quadratic equation a*x^2 + b*x + c = 0.
	D = b^2 - 4*a*c
	if (D >= 0) {	# Roots are real.
		x1 = (-b+sqrt(D))/(2*a)
		x2 = (-b-sqrt(D))/(2*a)
	}
	else {	# Roots are complex.
		x1 = complex(real=(-b/2/a), imaginary=(sqrt(-D)/2/a))
		x2 = complex(real=(-b/2/a), imaginary=(-sqrt(-D)/2/a))
	}
	return(c(x1, x2))
}
		
################################################################################

# Cross product of two vectors:
vec.cross = function(a, b) {
	# This function calculates the cross product of two vectors a x b.
	n = NULL
	n[1] = a[2]*b[3]-b[2]*a[3]
	n[2] = b[1]*a[3]-a[1]*b[3]
	n[3] = a[1]*b[2]-b[1]*a[2]
	return(n)
}

################################################################################

# Calculating centered values for a multidimensional array:
c.val = function(X) {
	# This function calculates the centered values (X - mean(X)) for data along the third dimension of a 3-dimensional array X.
	X.c = array(0, dim=dim(X))
	for (i in 1:dim(X)[1]) {
		for (j in 1:dim(X)[2]) {
			X.mean = mean(X[i,j,], na.rm=TRUE)
			X.c[i,j,] = X[i,j,]-X.mean
			}
		}
	return(X.c)
	}			

################################################################################

# Calculating the corresponding threshold R-value for a given p-value:
find.Rlim = function(pval, n) {
	# "pval" is the desired p-value; "n" is the sample size.
	t0 = qt(p=(1 - pval/2), df=(n - 1))  # The corresponding t-statistics
	R = sqrt(t0^2/(n - 2 + t0^2))        # The corresponding R-value (correlation coefficient).
	return(R)
}

################################################################################

# Test for multicollinearity using Variance Inflation Factor:
find.VIF = function(X) {
	# This function tests for the severity of the problem of multicollinearity based on VIF. VIF should not be > 5 for multicollinearity to be considered inconsequential.
	# "X" is a matrix with each column being an explanatory (predictor) variable.
	VIF = rep(0, times=ncol(X))
	for (i in 1:ncol(X)) {
		model = lm(X[,i] ~ X[,-i], na.action=na.exclude)
		adj.R2 = summary.lm(model)$adj.r.squared
		VIF[i] = 1/(1-adj.R2)
	}
	return(VIF)
}

################################################################################

# Multiple linear regression:
ml.reg = function(X, Y, plim=0.05, normalized=FALSE) {
	# This function performs multiple linear regression of Y on X, obtaining only statistically significant regression coefficients (p-value < "plim").
	# Y is a vector of response (dependent) variable and X is a matrix with columns of explanatory (independent) variables without the constant term.
	# "ml.reg" returns a matrix, first column being the regression coefficient (with the constant term as the first entry), second column being the standard error, and third column being the associated p-value.
	# If "normalized" is TRUE, then the inputs X and Y are normalized to their standard deviations.
	reg.coef = rep(0, times=(ncol(as.matrix(X))+1))
	se = rep(0, times=(ncol(as.matrix(X))+1))
	pval = rep(0, times=(ncol(as.matrix(X))+1))
	MAT = na.omit(cbind(Y, X))
	Y = MAT[,1]
	X = MAT[,-1]
	if (nrow(MAT) <= ncol(MAT)) {
		reg.coef = NaN
		se = NaN
		pval = NaN
	} else {
		if (normalized) {
			Y = (Y - mean(Y))/sd(Y)
			if (ncol(as.matrix(X)) == 1) X = (X - mean(X))/sd(X) else {
				for (j in 1:ncol(as.matrix(X))) X[,j] = (X[,j] - mean(X[,j]))/sd(X[,j])
			}
		}
		reg.model = lm(Y ~ X, na.action=na.exclude)
		reg.coef = summary(reg.model)$coefficients[,1]
		se = summary(reg.model)$coefficients[,2]
		pval = summary(reg.model)$coefficients[,4]
		pval[which(pval == "NaN")] = 1
		reg.coef[which(pval > plim)] = NaN
		se[which(pval > plim)] = NaN
	}
	if (length(reg.coef) != (ncol(as.matrix(X))+1)) {	# Singularities may have occured.
		reg.coef = rep(NaN, times=(ncol(as.matrix(X))+1))
		se = rep(NaN, times=(ncol(as.matrix(X))+1))
		pval = rep(NaN, times=(ncol(as.matrix(X))+1))
	}
	if (normalized) {
		reg.coef = reg.coef[-1]
		se = se[-1]
		pval = pval[-1]
	}
	return(cbind(reg.coef, se, pval))
}

################################################################################

# Orthogonal regression:
ortho.reg = function(X, Y, alpha=0.05) {
	# This function performs orthogonal regression of Y on X, i.e., minimizes the perpendicular distances between observations and regression line Y = b0 + b2*X.
	# X and Y are vectors of the explanatory (independent) and response (dependent) variables, respectively.
	# This function returns a vector c(b0, b1, b1.low, b1.up, R^2) for the linear model Y = b0 + b1*X, where b1.low and b1.up are the 100*(1-alpha)% confidence interval for b1; R^2 is the coefficient of determination.
	MAT = na.omit(cbind(X, Y))
	X = MAT[,1]; Y = MAT[,2]
	n = nrow(MAT)
	Sxx = sum((X-mean(X))^2); Syy = sum((Y-mean(Y))^2); Sxy = sum((X-mean(X))*(Y-mean(Y)));
	b1 = (Syy-Sxx+sqrt((Syy-Sxx)^2+4*Sxy^2))/2/Sxy
	b0 = mean(Y) - b1*mean(X)
	R = sum((X-mean(X))*(Y-mean(Y)))/sqrt(sum((X-mean(X))^2))/sqrt(sum((Y-mean(Y))^2))
	S = cbind(c(Sxx, Sxy), c(Sxy, Syy))
	LAMBDA = eigen(S)$values
	l1 = sort(LAMBDA)[2]; l2 = sort(LAMBDA)[1];
	theta = atan(b1)
	phi = asin((qchisq(1-alpha,1)/((n-1)*(l1/l2 + l2/l1 - 2)))^0.5)
	ci = c(tan(theta-phi), tan(theta+phi))
	return(c(b0, b1, ci, R^2))
}

################################################################################

# Multiple roots-finding:
multiroot = function(f, interval, tol=1e-9) {
	# This function searches for and calculates all possible roots of an arbitrary function within a given search interval.
	# "f" is the function to solve, e.g. f = function(x) x^2 - 3*x
	# "interval" is a vector of length 2, specifying the search interal.
	# "tol" is the tolerance level used in the function "uniroot".
	x = seq(interval[1], interval[2], length.out=1e3)
	y = rep(NaN, times=length(x))
	for (i in 1:length(y)) y[i] = f(x[i])
	count = 0
	root.int = NULL
	for (i in 1:(length(y)-1)) {
		# Case 1: no sign change
		if (sign(y[i+1]) == sign(y[i])) { } else {
		   # Case 2: sign change
			if (sign(y[i]) == 0) { } else {
				count = count + 1
				if (sign(y[i+1]) == 0) root.int[(count*2-1):(count*2)] = c(x[i], x[i+2]) else root.int[(count*2-1):(count*2)] = c(x[i], x[i+1])
			}
		}
	}
	roots = NULL
	for (i in 1:count) roots[i] = uniroot(f, root.int[(i*2-1):(i*2)], tol=tol)$root
	return(roots)
}
	
################################################################################
	
# Error functions:
# These are the Gauss error function, complementary error function, and the inverses.
# Error function (see Abramowitz and Stegun 29.2.29):
erf = function(x) 2*pnorm(x*sqrt(2)) - 1
# Complementary error function:
erfc = function(x) 2*pnorm(x*sqrt(2), lower=FALSE)
# The inverses of "erf" and "erfc":
erfinv = function(x) qnorm((1 + x)/2)/sqrt(2)
erfcinv = function(x) qnorm(x/2, lower=FALSE)/sqrt(2)

################################################################################

# Dicrete Fourier Series:
discrete.fourier = function(y, t) {
	# This function computes the discrete Fourier transformation for a time series y as a function of time domain t, assuming t is equally spaced and complete (i.e. no missing data). The function also returns the corresponding Fourier line spectrum in terms of R^2 = (n/2)*Ck^2/(n-1)/sy^2, where sy^2 is the variance of y.
	n = length(t)
	omega1 = 2*pi/n
	K = floor(n/2)
	A = rep(NaN, times=K)
	B = rep(NaN, times=K)
	C = rep(NaN, times=K)
	phi = rep(NaN, times=K)
	R.sq = rep(NaN, times=K)
	COS.MAT = matrix(NaN, nrow=n, ncol=K)
	SIN.MAT = matrix(NaN, nrow=n, ncol=K)
	for (k in 1:K) {
		if ((round(n/2)*2 == n) & k == K) {
			A[k] = 1/n*sum(y*cos(k*omega1*t))
			B[k] = 0/n*sum(y*sin(k*omega1*t))
		} else {
			A[k] = 2/n*sum(y*cos(k*omega1*t))
			B[k] = 2/n*sum(y*sin(k*omega1*t))
		}
		C[k] = sqrt(A[k]^2 + B[k]^2)
		if (A[k] > 0) phi[k] = atan(B[k]/A[k])
		if (A[k] < 0) phi[k] = atan(B[k]/A[k]) + pi
		if (A[k] == 0) phi[k] = pi/2
		R.sq[k] = n/2*C[k]^2/(n - 1)/var(y)
		COS.MAT[,k] = cos(k*omega1*t)
		SIN.MAT[,k] = sin(k*omega1*t)
	}
	y.fit = mean(y) + COS.MAT %*% A + SIN.MAT %*% B
	k = seq(1, K)
	omega = k*omega1
	fourier.fit = as.data.frame(cbind(k, omega, A, B, C, phi, R.sq))
	OUTPUT = list(n=n, y=y, var.y=var(y, na.rm=TRUE), y.fit=y.fit, cos.terms=COS.MAT, sin.terms=SIN.MAT, parameters=fourier.fit, R.sq=R.sq, C.sq=C^2)
	return(OUTPUT)
}

################################################################################

# Quantile of probability function y = f(x):
quantile.fx = function(x, y, p) {
	# This function finds the p-th quantile value of x for a given probability distribution vector y = f(x).
	# E.g. quantile.fx(x, y, 0.5) gives the median value of x given probability distribution y = f(x).
	cum.fx = cumsum(y)/tail(cumsum(y), 1)
	quant.x = spline(x=cum.fx, y=x, xout=p)$y
	return(quant.x)
}

################################################################################

# Mode of probability function y = f(x):
mode.fx = function(x, y) {
	# This function finds the mode of x given a probability distribution vector y = f(x), using the Golden Section method. In case of multi-modal distribution, the global maximum will be found.
	if (which(y == max(y)) == 1 | which(y == max(y)) == length(y)) {
		# Maximum is on the boundary.
		a = x[which(y == max(y))]
		print('Max is on boundary.')
	} else {
		tol = 1e-10
		tau = (sqrt(5) - 1)/2
		a = x[which(y == max(y))-1]
		b = x[which(y == max(y))+1]
		x1 = a + (1 - tau)*(b - a)
		y1 = spline(x, y, xout=x1)$y
		x2 = a + tau*(b - a)
		y2 = spline(x, y, xout=x2)$y
		while (abs(b - a) > tol) {
			if (y1 > y2) {
				a = x1
				x1 = x2
				y1 = y2
				x2 = a + tau*(b - a)
				y2 = spline(x, y, xout=x2)$y
			} else {
				b = x2
				x2 = x1
				y2 = y1
				x1 = a + (1 - tau)*(b - a)
				y1 = spline(x, y, xout=x1)$y
			}
		}
	}
	return(a)
}

################################################################################

# Find cumulative probability for where a specific value lies in a given sample distribution:
quantile.inv = function(x, x.spec) {
	# This function finds the cumulative probability for where a specific sample value "x.spec" lies within a given sample distribution "x", where "x" is a vector of sample values.
	# "x.spec" cannot be the exact median of "x".
	tol = 1e-9
	x = na.omit(x)
	if (x.spec > max(x) | x.spec < min(x)) {
		print('Requested sample value lies outside of the range of sample.')
		c = NaN
	} else {
		if (x.spec > median(x)) {
			a = 0.5
			b = 1
		} else {
			a = 0
			b = 0.5
		}
		c = a + (b - a)/2
		x1 = quantile(x, c)
		while (abs(x.spec - x1) > tol) {
			if (x.spec > x1) {
				a = c
				b = b
			} else {
				a = a
				b = c
			}
			c = a + (b - a)/2
			x1 = quantile(x, c)
		}
	}
	return(c)
}
			
################################################################################

# Trim data outside a specified quantile range:
trim.quant = function(X, quant=0.95, two.sided=TRUE) {
   # This function turns all data outside a specified quantile range into NA values.
   # For "two.sided=TRUE", "quant" specifies the quantile range to include. E.g., setting quant = 0.90 would trim data below the 0.05 quantile and above the 0.95 quantile, leaving only the "middle" 90% of data intact.
   # For "two.sided=FALSE", "quant" if < 0.5 specifies the quantile below which data are trimmed, and if > 0.5 specifies the quantile above which data are trimmed. E.g., setting quant = 0.95 would trim the top 5% of data.
   sub.data = X
   if (two.sided) {
      quant.range = quantile(X, probs=c((1 - quant)/2, (1 + quant)/2), na.rm=TRUE)
      sub.data[which(X <= quant.range[1] | X >= quant.range[2])] = NaN
   } else {
      quant.range = quantile(X, probs=quant, na.rm=TRUE)
      if (quant < 0.5) sub.data[which(X <= quant.range)] = NaN else sub.data[which(X >= quant.range)] = NaN
   }
   return(sub.data)
}

################################################################################

# Specify plot region:
# This function overrides default plot parameter.

quartz.par = function(title='', width=5.5, height=4, mfrow=c(1,1), mai=c(0.5,0.5,0.2,0.2), mgp=c(1.4,0.4,0), tcl=-0.25, ps=12) {
	quartz(title=title, width=width, height=height)
	par(mfrow=mfrow, mai=mai, mgp=mgp, tcl=tcl, ps=ps)
}

X11.par = function(title='', width=5.5, height=4, mfrow=c(1,1), mai=c(0.5,0.5,0.2,0.2), mgp=c(1.4,0.4,0), tcl=-0.25, ps=12) {
	X11(title=title, width=width, height=height)
	par(mfrow=mfrow, mai=mai, mgp=mgp, tcl=tcl, ps=ps)
}

################################################################################

# Standard major-axis regression:
sma.reg = function(X, Y, alpha=0.05) {
	# This function calculates the standard major-axis regression results for model Y = b0 + b1*X, giving also (1-alpha)% confidence intervals and r squared.
	MAT = na.omit(cbind(X, Y))
	X = MAT[,1]; Y = MAT[,2]
	n = length(X)
	r = cor(X, Y)
	b1 = sign(r)*sd(Y)/sd(X)
	b0 = mean(Y) - b1*mean(X)
	B = qt((1 - alpha/2), n-2)^2*(1 - r^2)/(n - 2)
	if (b1 > 0) {
		b1.low = b1*(sqrt(B + 1) - sqrt(B))
		b1.up = b1*(sqrt(B + 1) + sqrt(B))
	} else {
		b1.low = b1*(sqrt(B + 1) + sqrt(B))
		b1.up = b1*(sqrt(B + 1) - sqrt(B))
	}
	if (mean(X) >= 0) {
		b0.low = mean(Y) - b1.up*mean(X)
		b0.up = mean(Y) - b1.low*mean(X)
	} else {
		b0.low = mean(Y) - b1.low*mean(X)
		b0.up = mean(Y) - b1.up*mean(X)
	}
	out = list(reg.coef=c(b0, b1), conf.int=rbind(c(b0.low, b0.up), c(b1.low, b1.up)), r.squared=r^2)
	return(out)
}
		
################################################################################

# Model biases and errors:
model.bias.error = function(obs, mod) {
	# This function calculates various statistical measures used in model performance evaluation, given observed ('obs') and modeled ('mod') results.
	MAT = na.omit(cbind(obs, mod))
	obs = MAT[,1]; mod = MAT[,2]
	n = length(obs)
	# Mean bias:
	MB = sum(mod - obs)/n
	# Normalized mean bias:
	NMB = sum(mod - obs)/sum(obs)*100
	# Normalized mean error:
	NME = sum(abs(mod - obs))/sum(obs)*100
	# Root mean square error:
	RMSE = sqrt(sum((mod - obs)^2)/n)
	out = list(MB=MB, NMB=NMB, NME=NME, RMSE=RMSE)
	return(out)
}

################################################################################

# Scatterplot of two variables with statistics:

plot.reg = function(x, y, x.quant=1, y.quant=1, two.sided=TRUE, top.left=TRUE, digits=3, xlab='x', ylab='y', mai=c(0.4, 0.4, 0.2, 0.2), mgp=c(1.5, 0.5, 0), tcl=-0.25, ps=14, cex=1.1) {
   # This function plots a scatterplot of two variables, x and y, also displaying their correlation coefficient and regression results.
   # It has the functionality to trim the data first according to a specified interquantile range to remain, using the function "trim.quant()".
   # "x.quant" and "y.quant" specify the quantile range of x and y to remain, outside of which the data will be trimmed. See "trim.quant()" for the details.
   # If "topleft=TRUE", the statistics will be printed on the top left corner of the scatterplot. Otherwise, they will be on the bottom left corner.
   
   # Trimming data:
   if (x.quant < 1)  x = trim.quant(x, quant=x.quant, two.sided=two.sided)
   if (y.quant < 1)  y = trim.quant(y, quant=y.quant, two.sided=two.sided)
   
   par(mai=mai, mgp=mgp, tcl=tcl, ps=ps)
   
   x = as.vector(x)
   y = as.vector(y)
   plot(x, y, pch=21, xlab=xlab, ylab=ylab)
   cor.xy = cor(x, y, use="complete.obs")
   reg = lm(y ~ x)
   abline(reg, col="blue", lty=2)
   b0 = summary(reg)$coefficients[1,1]
   b1 = summary(reg)$coefficients[2,1]
   pval = summary(reg)$coefficients[2,4]
   
   xlim = range(x, na.rm=TRUE)
   ylim = range(y, na.rm=TRUE)
   xpos = xlim[1] + 0.010*diff(xlim)
   if (top.left) {
      ypos1 = ylim[1] + 0.980*diff(ylim)
      ypos2 = ylim[1] + 0.890*diff(ylim)
      ypos3 = ylim[1] + 0.800*diff(ylim)
   } else {
      ypos1 = ylim[1] + 0.230*diff(ylim)
      ypos2 = ylim[1] + 0.140*diff(ylim)
      ypos3 = ylim[1] + 0.050*diff(ylim)
   }
   text(x=xpos, y=ypos1, labels=paste('r  = ', as.character(signif(cor.xy, digits=digits)), sep=''), cex=cex, col='blue', adj=c(0,1))
   text(x=xpos, y=ypos2, labels=paste('y = ', as.character(signif(b0, digits=digits)), ' + ', as.character(signif(b1, digits=digits)), ' x', sep=''), cex=cex, col='blue', adj=c(0,1))
   if (pval < 1e-99) {
      text(x=xpos, y=ypos3, labels=paste('p < ', as.character(1e-99), sep=''), cex=cex, col='blue', adj=c(0,1))
   } else {
      text(x=xpos, y=ypos3, labels=paste('p = ', as.character(signif(pval, digits=digits)), sep=''), cex=cex, col='blue', adj=c(0,1))
   }
   
   return(reg)
   
}

################################################################################

# Principal component analysis:
find.PC = function(X, ij=NULL, normalized=TRUE, date.vec=NULL, var.name=NULL) {
   # This function performs PCA on a given data array X.
   # "ij" specifies what type (dimension) of data array X is. If ij = NULL, X is simply a usual 2D data matrix where rows = observations (n; usually times) and columns = variables (K). If ij = c(i,j), X is a 4D array where the first two dimensions are spatial coordinates (e.g., dim1 = lon, dim2 = lat); then i and j are indexing the exact location for which the PCA is to be done, and the remaining two dimensions carry the same meaning as a usual 2D data matrix. If ij = rbind(c(i1,j1), c(i2,j2), ...), X is still a 4D data array with the same meaning, but we are interested in more than one locations in the array, and the averages over all the locations are to be computed before PCA is done.
   # "normalized" sets whether the data should be centered and standardized before PCA. This is highly recommended if the variables have different units.
   # "date.vec" is a numeric vector of dates or times referring to the n observations.
   # "var.name" is a character vector of names of the K variables.
   if (is.null(ij)) {
      X.data = X[,]
   } else if (length(ij) == 2) {
      i = ij[1]; j = ij[2]
      X.data = X[i,j,,]
   } else {
      X.data = matrix(NaN, nrow=dim(X)[3], ncol=dim(X)[4])
      for (n in 1:nrow(X.data)) {
         for (k in 1:ncol(X.data)) {
            X.data[n,k] = mean(X[,,n,k][ij], na.rm=TRUE)
         }
      }
   }
   # Means and standard deviations of original data:
   X.mean = apply(X.data, MARGIN=2, FUN=mean, na.rm=TRUE)
   X.sd = apply(X.data, MARGIN=2, FUN=sd, na.rm=TRUE)
   # Find PC:
   pca = prcomp(~X.data, scale=normalized, na.action=na.exclude)
   PC.sd = pca$sdev
   PC.rot = pca$rotation
   PC.data = pca$x
   # Output:
   OUTPUT = list(date=date.vec, X.mean=X.mean, X.sd=X.sd, var.name=var.name, PC.sd=PC.sd, PC.rot=PC.rot, PC.data=PC.data)
   return(OUTPUT)
}

################################################################################

# Multiple linear regression with forward selection:
mlr.forward = function(y, X, dR2.lim=0.05, VIF.lim=5) {
   # This function performs multiple linear regression with forward selection based on the increases in adjusted R^2 value, and returns the indices for the predictor variables to be included in the regression model in their order of inclusion.
   # "y" is a vector of predictand.
   # "X" is a 2D data matrix where rows = observations (n) and columns = predictor variables (K).
   # "dR2.lim" sets the threshold of increased adjusted R^2 below which new variables will not be added to the model.
   # "VIF.lim" sets the threshold of variance inflation factor above which a new variable will not be added to the model. It utilizes another function "find.VIF" in the same file "get_stat.R".
   K_all = 1:ncol(X)
   cor_1st = NULL
   for (k in K_all) cor_1st = c(cor_1st, cor(y, X[,k], use='complete.obs'))
   K_inc = which(abs(cor_1st) == max(abs(cor_1st), na.rm=TRUE))
   K_next = K_all[-K_inc]
   repeat {
      dR2 = NULL
      for (k in K_next) {
         # Exclude variable that leads to VIF > VIF.lim:
         VIF = max(find.VIF(X[,c(K_inc, k)]), na.rm=TRUE)
         if (VIF > VIF.lim) dR2 = c(dR2, 0) else {
            mlr1 = lm(y ~ X[,c(K_inc)])
            mlr2 = lm(y ~ X[,c(K_inc, k)])
            dR2 = c(dR2, summary(mlr2)$adj.r.squared - summary(mlr1)$adj.r.squared)
         }
      }
      if (length(dR2) == 0) break else {
         dR2_max = max(dR2, na.rm=TRUE)
         if (dR2_max < dR2.lim) break else {
            K_inc = c(K_inc, K_next[which(dR2 == dR2_max)])
            K_next = K_all[-K_inc]
         }
      }
   }
   print('Variables included: ', quote=FALSE)
   print(K_inc)
   return(K_inc)
}

################################################################################
