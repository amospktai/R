##################################################################
##################################################################



# Plot NA vs. 2x2.5 and 4x5 regression coefficients:

# Model-Observation comparison:

plot.reg.compare.mapx4 = function(species, met, step=FALSE, std=FALSE, normalized=FALSE, norm.method='sd', include.plot=TRUE, zmax=NULL, mai=c(0.1, 0.1, 0.2, 0.1), tcl=-0.1, ps=12, legend.mar=4.8, legend.width=1.2, line=0.1, cex=1.2) {
	
	# Meteorological variables:
	if (met == "T") k = 1
	if (met == "RH") k = 2
	if (met == "prec") k = 3
	if (met == "hght") k = 4
	if (met == "dp/dt") k = 5
	if (met == "wind") k = 6
	if (met == "EW") k = 7
	if (met == "NS") k = 8

	# Species:
	if (species == "PM25") spec.GC = "PMtotchm"
	if (species == "SO4") spec.GC = "sulfate"
	if (species == "NO3") spec.GC = "nitrate"
	if (species == "OC") spec.GC = "OC"
	if (species == "EC") spec.GC = "BC"
	if (species == "SOA") spec.GC = "SOA"
	if (species == "POA") spec.GC = "POA"
	
	if (step == TRUE) {
		reg.mod = "step1"
		reg.obs = "step"
	} else if (std == TRUE) {
		reg.mod = "std"
		reg.obs = "std"
	} else {
		reg.mod = "fix"
		reg.obs = "fix"
	}

	# Modeled correlations:
	
	# NA nested grid:
	filename = paste("~/Dropbox/Research/PM_geos5/data_reg_NA_8x/reg_", reg.mod, "_", spec.GC, ".RData", sep="")
	data1 = load.RData(filename, "reg.out")[,,k]
	if (k == 5) data1 = data1/1000
	if (normalized) {
		model = 'NA'
		if (std) stop('Standardized data need not be normalized to mean.')
		if (norm.method == 'mean') {
			start.date = year*10000 + 101
			end.date = year*10000 + 1231
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/dim.RData", sep="")
			date.conc = load.RData(filename, 'date.vec')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/conc_", spec.GC, ".RData", sep="")
			conc = load.RData(filename, spec.GC)[,,which(date.conc >= start.date & date.conc <= end.date)]
			conc = apply(conc, c(1,2), mean, na.rm=TRUE)
			data1 = data1/conc*100
		}
		if (norm.method == 'sd') {
			specname = paste(spec.GC, '.std', sep='')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/conc_", spec.GC, ".RData", sep="")
			conc = load.RData(filename, specname)
			conc.sd = apply(conc, c(1,2), sd, na.rm=TRUE)
			data1 = data1/conc.sd
		}
	}		

	# 2x2.5:
	filename = paste("~/Dropbox/Research/PM_geos5/data_reg_2x2.5_201011/reg_", reg.mod, "_", spec.GC, ".RData", sep="")
	data2 = load.RData(filename, "reg.out")[,,k]
	if (k == 5) data2 = data2/1000
	if (normalized) {
		model = '2x2.5_201011'
		if (std) stop('Standardized data need not be normalized to mean.')
		if (norm.method == 'mean') {
			start.date = year*10000 + 101
			end.date = year*10000 + 1231
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/dim.RData", sep="")
			date.conc = load.RData(filename, 'date.vec')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/conc_", spec.GC, ".RData", sep="")
			conc = load.RData(filename, spec.GC)[,,which(date.conc >= start.date & date.conc <= end.date)]
			conc = apply(conc, c(1,2), mean, na.rm=TRUE)
			data2 = data2/conc*100
		}
		if (norm.method == 'sd') {
			specname = paste(spec.GC, '.std', sep='')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/conc_", spec.GC, ".RData", sep="")
			conc = load.RData(filename, specname)
			conc.sd = apply(conc, c(1,2), sd, na.rm=TRUE)
			data2 = data2/conc.sd
		}
	}		
	
	# 4x5:
	filename = paste("~/Dropbox/Research/PM_geos5/data_reg_4x5_201011/reg_", reg.mod, "_", spec.GC, ".RData", sep="")
	data3 = load.RData(filename, "reg.out")[,,k]
	if (k == 5) data3 = data3/1000
	if (normalized) {
		model = '4x5_201011'
		if (std) stop('Standardized data need not be normalized to mean.')
		if (norm.method == 'mean') {
			start.date = year*10000 + 101
			end.date = year*10000 + 1231
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/dim.RData", sep="")
			date.conc = load.RData(filename, 'date.vec')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/conc_", spec.GC, ".RData", sep="")
			conc = load.RData(filename, spec.GC)[,,which(date.conc >= start.date & date.conc <= end.date)]
			conc = apply(conc, c(1,2), mean, na.rm=TRUE)
			data3 = data3/conc*100
		}
		if (norm.method == 'sd') {
			specname = paste(spec.GC, '.std', sep='')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/conc_", spec.GC, ".RData", sep="")
			conc = load.RData(filename, specname)
			conc.sd = apply(conc, c(1,2), sd, na.rm=TRUE)
			data3 = data3/conc.sd
		}
	}		

	
	# Observed correlations:
	
	# EPA-AQS / GEOS5:
	specname = species
	if (species == "SOA" | species == "POA") specname = "OC"
	filename = paste("~/Dropbox/Research/PM_geos5/data_reg_AQS_GEOS5/reg_", reg.obs, "_", specname, ".RData", sep="")
	data4 = rm.nonUS.2x2.5(load.RData(filename, "reg.out")[,,k])
	if (k == 5) data4 = data4/1000
	if (normalized) {
		if (std) stop('Standardized data need not be normalized.')
		if (norm.method == 'mean') {
			spec.AQS = paste(species, '.invdist', sep='')
			start.date = year*10000 + 101
			end.date = year*10000 + 1231
			filename = '~/Dropbox/Research/EPA_data/data_conc_AQS/AQS_GRID_2x2.5_2004-2008.RData'
			date.conc = load.RData(filename, 'date.vec')
			conc = load.RData(filename, spec.AQS)[,,which(date.conc >= start.date & date.conc <= end.date)]
			conc = apply(conc, c(1,2), mean, na.rm=TRUE)
			data4 = data4/conc*100
		}
		if (norm.method == 'sd') {
			spec.AQS = paste(species, '.std', sep='')
			filename = '~/Dropbox/Research/EPA_data/data_conc_AQS/AQS_GRID_2x2.5_2004-2008.RData'
			conc = load.RData(filename, spec.AQS)
			conc.sd = apply(conc, c(1,2), sd, na.rm=TRUE)
			data4 = data4/conc.sd
		}
	}		
	
	if (include.plot) {
		
		# Decide zmax:
		if (is.null(zmax)) zmax = quantile(c(abs(as.vector(data1)), abs(as.vector(data2)), abs(as.vector(data3)), abs(as.vector(data4))), probs=0.975, na.rm=TRUE)
		data1[which(abs(data1) > zmax)] = sign(data1[which(abs(data1) > zmax)])*zmax
		data2[which(abs(data2) > zmax)] = sign(data2[which(abs(data2) > zmax)])*zmax
		data3[which(abs(data3) > zmax)] = sign(data3[which(abs(data3) > zmax)])*zmax
		data4[which(abs(data4) > zmax)] = sign(data4[which(abs(data4) > zmax)])*zmax
	
		# Plotting:	
		par(mfrow=c(2,2))

		filename = "~/Dropbox/Research/PM_geos5/data_conc_NA/dim.RData"
		lon.plot = load.RData(filename, "lon.US")
		lat.plot = load.RData(filename, "lat.US")
		plot.field(data1, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
		add.text = expression('GEOS-Chem '*0.5*degree*'x'*0.667*degree)
		mtext(add.text, side=3, line=line, cex=cex)
	
		filename = "~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/dim.RData"
		lon.plot = load.RData(filename, "lon.US")
		lat.plot = load.RData(filename, "lat.US")
		plot.field(data2, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
		add.text = expression('GEOS-Chem '*2*degree*'x'*2.5*degree)
		mtext(add.text, side=3, line=line, cex=cex)
	
		filename = "~/Dropbox/Research/PM_geos5/data_conc_4x5_201011/dim.RData"
		lon.plot = load.RData(filename, "lon.US")
		lat.plot = load.RData(filename, "lat.US")
		plot.field(data3, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
		add.text = add.text = expression('GEOS-Chem '*4*degree*'x'*5*degree)
		mtext(add.text, side=3, line=line, cex=cex)
	
		filename = "~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/dim.RData"
		lon.plot = load.RData(filename, "lon.US")
		lat.plot = load.RData(filename, "lat.US")
		plot.field(data4, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
		add.text = "EPA-AQS observations"
		mtext(add.text, side=3, line=line, cex=cex)
		
	}
	
	OUTPUT = list(mod.NA=data1, mod.2x2.5=data2, mod.4x5=data3, obs=data4)
	return(OUTPUT)
	
}



##################################################################



# Compare the regression coefficients across model resolutions:

plot.reg.compare.scat = function(species, met, resol1, resol2=NULL, obs=FALSE, step=FALSE, std=FALSE, normalized=FALSE, norm.method='sd', cex=1.6) {
	
	reg.data = plot.reg.compare.mapx4(species, met, step=step, std=std, normalized=normalized, norm.method=norm.method, include.plot=FALSE)
	
	if (!is.null(resol2)) {
		
		# Modeled correlations:
		
		if (resol1 == 'NA') {
			data1 = reg.data$mod.NA
			filename = "~/Dropbox/Research/PM_geos5/data_conc_NA/dim.RData"
			lon1 = load.RData(filename, "lon.US")
			lat1 = load.RData(filename, "lat.US")
		}
		if (resol1 == '2x2.5') {
			data1 = rm.nonUS.2x2.5(reg.data$mod.2x2.5)
			filename = "~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/dim.RData"
			lon1 = load.RData(filename, "lon.US")
			lat1 = load.RData(filename, "lat.US")
		}
		if (resol2 == '2x2.5') {
			data2 = rm.nonUS.2x2.5(reg.data$mod.2x2.5)
			filename = "~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/dim.RData"
			lon2 = load.RData(filename, "lon.US")
			lat2 = load.RData(filename, "lat.US")
		}
		if (resol2 == '4x5') {
			data2 = reg.data$mod.4x5
			filename = "~/Dropbox/Research/PM_geos5/data_conc_4x5_201011/dim.RData"
			lon2 = load.RData(filename, "lon.US")
			lat2 = load.RData(filename, "lat.US")
		}
				
		X = as.vector(sp.dissolve(data1[,length(lat1):1], lon1, rev(lat1), lon2, rev(lat2))[,length(lat2):1])
		Y = as.vector(data2)
			
		reg = lmodel2(Y ~ X)
		r.sq = reg$rsquare
		a = reg$regression.results[2,2]
		b = reg$regression.results[2,3]
		xylim = c(min(c(X, Y), na.rm=TRUE), max(c(X, Y), na.rm=TRUE))
		plot(X, Y, xlim=xylim, ylim=xylim, xlab='', ylab='')
		abline(a=0, b=1)
		abline(a=a, b=b, lty=2)
		text(x=(xylim[1] + 0.060*diff(xylim)), y=(xylim[1] + 0.961*diff(xylim)), labels=expression(R^2), cex=cex)
		text(x=(xylim[1] + 0.255*diff(xylim)), y=(xylim[1] + 0.950*diff(xylim)), labels=paste('= ', as.character(round(r.sq, digits=1.5)), sep=''), cex=cex)
		text(x=(xylim[1] + 0.220*diff(xylim)), y=(xylim[1] + 0.850*diff(xylim)), labels=paste('b = ', as.character(round(b, digits=3)), sep=''), cex=cex)
		#text(x=(xylim[1] + 0.85*diff(xylim)), y=(xylim[1] + 0.15*diff(xylim)), labels=species, cex=2.0)
		
	}
		
	if (is.null(resol2) & obs) {
		
		# Only allow comparison for 2x2.5 resolution.
		if (resol1 != '2x2.5') stop('Model resolution specified is not allowed!')
			
		# Modeled correlations:
		# 2x2.5:
		data1 = rm.nonUS.2x2.5(reg.data$mod.2x2.5)
		
		# Observed correlations:
		# EPA-AQS / GEOS5:
		data2 = rm.nonUS.2x2.5(reg.data$obs)
	
		X = as.vector(data2)
		Y = as.vector(data1)
			
		reg = lmodel2(Y ~ X)
		r.sq = reg$rsquare
		a = reg$regression.results[2,2]
		b = reg$regression.results[2,3]
		xylim = c(min(c(X, Y), na.rm=TRUE), max(c(X, Y), na.rm=TRUE))
		plot(X, Y, xlim=xylim, ylim=xylim, xlab='', ylab='')
		abline(a=0, b=1)
		abline(a=a, b=b, lty=2)
		text(x=(xylim[1] + 0.060*diff(xylim)), y=(xylim[1] + 0.961*diff(xylim)), labels=expression(R^2), cex=cex)
		text(x=(xylim[1] + 0.255*diff(xylim)), y=(xylim[1] + 0.950*diff(xylim)), labels=paste('= ', as.character(round(r.sq, digits=1.5)), sep=''), cex=cex)
		text(x=(xylim[1] + 0.220*diff(xylim)), y=(xylim[1] + 0.850*diff(xylim)), labels=paste('b = ', as.character(round(b, digits=3)), sep=''), cex=cex)
		#text(x=(xylim[1] + 0.85*diff(xylim)), y=(xylim[1] + 0.15*diff(xylim)), labels=species, cex=2.0)

	}
	
	return(reg)
	
}



##################################################################



# Comparing 2x2.5 simulated vs. observed regression coefficients:

plot.reg.compare.mapx2 = function(species, met, model="2x2.5_201011", step=FALSE, std=FALSE, normalized=FALSE, norm.method='sd', year=2006, include.plot=TRUE, zmax=NULL, mai=c(0.2, 0.2, 0.2, 0.2), tcl=-0.2, ps=12, legend.mar=4.8, legend.width=1.5, cex=1.2, line=0.1, mtext1='Simulated reg coef', mtext2='Observed reg coef', season=NULL, model.left=TRUE) {
	
	# Meteorological variables:
	if (met == "T") { k = 1; met.name = 'TMPU' }
	if (met == "RH") { k = 2; met.name = 'RH' }
	if (met == "prec") { k = 3; met.name = 'PREACC' }
	if (met == "hght") { k = 4; met.name = 'Z850' }
	if (met == "dp/dt") { k = 5; met.name = 'DSLPDT' }
	if (met == "wind") { k = 6; met.name = 'WINDS' }
	if (met == "EW") { k = 7; met.name = 'EWDIR' }
	if (met == "NS") { k = 8; met.name = 'NSDIR' }

	# Species:
	if (species == "PM25") spec.GC = "PMtotchm"
	if (species == "SO4") spec.GC = "sulfate"
	if (species == "NO3") spec.GC = "nitrate"
	if (species == "OC") spec.GC = "OC"
	if (species == "EC") spec.GC = "BC"
	if (species == "SOA") spec.GC = "SOA"
	if (species == "POA") spec.GC = "POA"
	
	if (step == TRUE) {
		reg.mod = "step1"
		reg.obs = "step"
	} else if (std == TRUE) {
		reg.mod = "std"
		reg.obs = "std"
	} else {
		reg.mod = "fix"
		reg.obs = "fix"
	}
	
	# Modeled correlations:
	# Determined by "model":
	if (is.null(season)) filename = paste("~/Dropbox/Research/PM_geos5/data_reg_", model, "/reg_", reg.mod, "_", spec.GC, ".RData", sep="")
	else filename = paste("~/Dropbox/Research/PM_geos5/data_reg_", model, "/reg_", reg.mod, "_", season, "_", spec.GC, ".RData", sep="")
	data1 = rm.nonUS.2x2.5(load.RData(filename, "reg.out")[,,k])
	if (k == 5) data1 = data1/1000
	
	if (normalized) {
		if (std) stop('Standardized data need not be normalized to mean.')
		if (norm.method == 'mean') {
			start.date = year*10000 + 101
			end.date = year*10000 + 1231
			if (!is.null(season)) stop('Seasonal subsets are not allowed for displaying normalized results.')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/dim.RData", sep="")
			date.conc = load.RData(filename, 'date.vec')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/conc_", spec.GC, ".RData", sep="")
			conc = load.RData(filename, spec.GC)[,,which(date.conc >= start.date & date.conc <= end.date)]
			conc = apply(conc, c(1,2), mean, na.rm=TRUE)
			data1 = data1/conc*100
		}
		if (norm.method == 'sd') {
			if (!is.null(season)) stop('Seasonal subsets are not allowed for displaying normalized results.')
			specname = paste(spec.GC, '.std', sep='')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/conc_", spec.GC, ".RData", sep="")
			conc = load.RData(filename, specname)
			conc.sd = apply(conc, c(1,2), sd, na.rm=TRUE)
			data1 = data1/conc.sd
		}
	}		
	
	# Observed correlations:
	# EPA-AQS / GEOS5:
	specname = species
	if (species == "SOA" | species == "POA") specname = "OC"
	if (is.null(season)) filename = paste("~/Dropbox/Research/PM_geos5/data_reg_AQS_GEOS5/reg_", reg.obs, "_", specname, ".RData", sep="")
	else filename = paste("~/Dropbox/Research/PM_geos5/data_reg_AQS_GEOS5/reg_", reg.obs, "_", season, "_", specname, ".RData", sep="")
	data2 = rm.nonUS.2x2.5(load.RData(filename, "reg.out")[,,k])
	if (k == 5) data2 = data2/1000
	
	if (normalized) {
		if (std) stop('Standardized data need not be normalized.')
		if (norm.method == 'mean') {
			spec.AQS = paste(species, '.invdist', sep='')
			start.date = year*10000 + 101
			end.date = year*10000 + 1231
			if (!is.null(season)) stop('Seasonal subsets are not allowed for displaying normalized results.')
			filename = '~/Dropbox/Research/EPA_data/data_conc_AQS/AQS_GRID_2x2.5_2004-2008.RData'
			date.conc = load.RData(filename, 'date.vec')
			conc = load.RData(filename, spec.AQS)[,,which(date.conc >= start.date & date.conc <= end.date)]
			conc = apply(conc, c(1,2), mean, na.rm=TRUE)
			data2 = data2/conc*100
		}
		if (norm.method == 'sd') {
			spec.AQS = paste(species, '.std', sep='')
			if (!is.null(season)) stop('Seasonal subsets are not allowed for displaying normalized results.')
			filename = '~/Dropbox/Research/EPA_data/data_conc_AQS/AQS_GRID_2x2.5_2004-2008.RData'
			conc = load.RData(filename, spec.AQS)
			conc.sd = apply(conc, c(1,2), sd, na.rm=TRUE)
			data2 = data2/conc.sd
		}
	}		
	
	if (include.plot) {
		
		# Decide zmax:
		if (is.null(zmax)) zmax = quantile(c(abs(as.vector(data1)), abs(as.vector(data2))), probs=0.975, na.rm=TRUE)
		data1[which(abs(data1) > zmax)] = sign(data1[which(abs(data1) > zmax)])*zmax
		data2[which(abs(data2) > zmax)] = sign(data2[which(abs(data2) > zmax)])*zmax

		# Plotting:
	
		filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/dim.RData", sep="")
		lon.plot = load.RData(filename, "lon.US")
		lat.plot = load.RData(filename, "lat.US")
		
		if (model.left) {
			
			plot.field(data1, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
			mtext(mtext1, cex=cex, line=line)

			plot.field(data2, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
			mtext(mtext2, cex=cex, line=line)
			
		} else {
			
			plot.field(data2, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
			mtext(mtext2, cex=cex, line=line)
			
			plot.field(data1, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
			mtext(mtext1, cex=cex, line=line)

		}
			
		
	}
	
	OUTPUT = list(mod=data1, obs=data2)
	return(OUTPUT)
	
}



##################################################################



# Plotting simulated correlations only:

plot.reg.compare.mapx1 = function(species, met, model="2x2.5_201011", region='US', China.only=FALSE, step=FALSE, std=FALSE, normalized=FALSE, norm.method='sd', year=2006, include.plot=TRUE, zmax=NULL, mai=c(0.2, 0.2, 0.2, 0.2), tcl=-0.2, ps=12, legend.mar=4.8, legend.width=1.5, cex=1.2, line=0.1, mtext1='Simulated reg coef', mtext2='Observed reg coef', season=NULL) {
	
	# Meteorological variables:
	if (met == "T") { k = 1; met.name = 'TMPU' }
	if (met == "RH") { k = 2; met.name = 'RH' }
	if (met == "prec") { k = 3; met.name = 'PREACC' }
	if (met == "hght") { k = 4; met.name = 'Z850' }
	if (met == "dp/dt") { k = 5; met.name = 'DSLPDT' }
	if (met == "wind") { k = 6; met.name = 'WINDS' }
	if (met == "EW") { k = 7; met.name = 'EWDIR' }
	if (met == "NS") { k = 8; met.name = 'NSDIR' }

	# Species:
	if (species == "PM25") spec.GC = "PMtotchm"
	if (species == "SO4") spec.GC = "sulfate"
	if (species == "NO3") spec.GC = "nitrate"
	if (species == "OC") spec.GC = "OC"
	if (species == "EC") spec.GC = "BC"
	if (species == "SOA") spec.GC = "SOA"
	if (species == "POA") spec.GC = "POA"
	
	if (step == TRUE) {
		reg.mod = "step1"
		reg.obs = "step"
	} else if (std == TRUE) {
		reg.mod = "std"
		reg.obs = "std"
	} else {
		reg.mod = "fix"
		reg.obs = "fix"
	}
	
	if (region == 'US') model = model
	if (region == 'EA') model = paste(model, '_EA', sep='')
	
	# Modeled correlations:
	# Determined by "model":
	if (is.null(season)) filename = paste("~/Dropbox/Research/PM_geos5/data_reg_", model, "/reg_", reg.mod, "_", spec.GC, ".RData", sep="")
	else filename = paste("~/Dropbox/Research/PM_geos5/data_reg_", model, "/reg_", reg.mod, "_", season, "_", spec.GC, ".RData", sep="")
	data1 = load.RData(filename, "reg.out")[,,k]
	if (region == 'US') data1 = rm.nonUS.2x2.5(data1)
	if (k == 5) data1 = data1/1000
	
	if (normalized) {
		if (std) stop('Standardized data need not be normalized to mean.')
		if (norm.method == 'mean') {
			start.date = year*10000 + 101
			end.date = year*10000 + 1231
			if (!is.null(season)) stop('Seasonal subsets are not allowed for displaying normalized results.')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/dim.RData", sep="")
			date.conc = load.RData(filename, 'date.vec')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/conc_", spec.GC, ".RData", sep="")
			conc = load.RData(filename, spec.GC)[,,which(date.conc >= start.date & date.conc <= end.date)]
			conc = apply(conc, c(1,2), mean, na.rm=TRUE)
			data1 = data1/conc*100
		}
		if (norm.method == 'sd') {
			if (!is.null(season)) stop('Seasonal subsets are not allowed for displaying normalized results.')
			specname = paste(spec.GC, '.std', sep='')
			filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/conc_", spec.GC, ".RData", sep="")
			conc = load.RData(filename, specname)
			conc.sd = apply(conc, c(1,2), sd, na.rm=TRUE)
			data1 = data1/conc.sd
		}
	}		
	
	if (China.only) {
		if (region != 'EA') stop('Cannot plot for China when region specified is not East Asia (EA).')
		filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/dim.RData", sep="")
		lon.plot = load.RData(filename, paste('lon.', region, sep=''))
		lat.plot = load.RData(filename, paste('lat.', region, sep=''))
		data1 = data1[which(lon.plot >= 72.5 & lon.plot <= 135),which(lat.plot >= 18 & lat.plot <= 54)]
	}
	
	if (include.plot) {
		
		# Decide zmax:
		if (is.null(zmax)) zmax = quantile(abs(as.vector(data1)), probs=0.975, na.rm=TRUE)
		data1[which(abs(data1) > zmax)] = sign(data1[which(abs(data1) > zmax)])*zmax
		
		# Plotting:
	
		filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", model, "/dim.RData", sep="")
		lon.plot = load.RData(filename, paste('lon.', region, sep=''))
		lat.plot = load.RData(filename, paste('lat.', region, sep=''))
		if (China.only) {
			lon.plot = lon.plot[which(lon.plot >= 72.5 & lon.plot <= 135)]
			lat.plot = lat.plot[which(lat.plot >= 18 & lat.plot <= 54)]
		}
		plot.field(data1, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
		mtext(mtext1, cex=cex, line=line)
		
	}
	
	return(data1)
	
}



##################################################################
##################################################################



# 6-panel comparison with:
# 1. GEOS-Chem NA nested grid
# 2. GEOS-Chem 2x3
# 3. GEOS-Chem 4x5
# 4. AQS-GEOS5
# 5. AQS-NCEP/NCAR
# 6. IMPROVE-NCEP/NCAR


# Model-Observation comparison:

plot.reg.compare.mapx6 = function(species, met, step=TRUE, normalized=FALSE, year=2006, zmax=NULL) {
	
	setwd('~/Dropbox/Research')
	
	# Meteorological variables:

	if (met == "T") k = 1
	if (met == "RH") k = 2
	if (met == "precip") k = 3
	if (met == "height") k = 4
	if (met == "dSLP/dt") k = 5
	if (met == "wind") k = 6
	if (met == "EWdir") k = 7
	if (met == "NSdir") k = 8

	# Species:

	if (species == "PM25") spec.GC = "PMtotchm"
	if (species == "SO4") spec.GC = "sulfate"
	if (species == "NO3") spec.GC = "nitrate"
	if (species == "OC") spec.GC = "OC"
	if (species == "EC") spec.GC = "BC"
	if (species == "SOA") spec.GC = "SOA"
	if (species == "POA") spec.GC = "POA"
	
	if (step == TRUE) {
		reg.mod = "step1"
		reg.obs = "step"
	}
	else {
		reg.mod = "fix"
		reg.obs = "fix"
	}

	# Modeled correlations:
	
	# NA nested grid:
	filename = paste("PM_geos5/data_reg_NA_8x/reg_", reg.mod, "_", spec.GC, ".RData", sep="")
	data1 = load.RData(filename, "reg.out")[,,k]
	if (k == 5) data1 = data1/1000
	
	# 2x2.5:
	filename = paste("PM_geos5/data_reg_2x2.5_201011/reg_", reg.mod, "_", spec.GC, ".RData", sep="")
	data2 = load.RData(filename, "reg.out")[,,k]
	if (k == 5) data2 = data2/1000
	
	# 4x5:
	filename = paste("PM_geos5/data_reg_4x5_201011/reg_", reg.mod, "_", spec.GC, ".RData", sep="")
	data3 = load.RData(filename, "reg.out")[,,k]
	if (k == 5) data3 = data3/1000
	
	# Observed correlations:
	
	# EPA-AQS / GEOS5:
	specname = species
	if (species == "SOA" | species == "POA") specname = "OC"
	filename = paste("PM_geos5/data_reg_AQS_GEOS5/reg_", reg.obs, "_", specname, ".RData", sep="")
	data4 = rm.nonUS.2x2.5(load.RData(filename, "reg.out")[,,k])
	if (k == 5) data4 = data4/1000
	
	# EPA-AQS / NCEP:
	specname = species
	if (species == "SOA" | species == "POA") specname = "OC"
	filename = paste("EPA_data/data_reg_AQS_NCEP/reg_step_", specname, ".RData", sep="")
	if (species == "PM25") specname = "invdist"
	dataname = paste("PM25.", specname, ".met.reg.step13", sep="")
	data5 = rm.nonUS(load.RData(filename, dataname)[,,k])
	data5 = data5[,dim(data5)[2]:1]
	if (k == 3) data5 = data5/25.4
	if (k == 6) data5 = data5*100
		
	# IMPROVE / NCEP:
	specname = species
	if (species == "SOA" | species == "POA") specname = "OC"
	if (species == "PM25") specname = spec.GC
	filename = paste("EPA_data/data_reg_IMPROVE_NCEP/reg_step_", specname, ".RData", sep="")
	data6 = rm.nonUS(load.RData(filename, "reg.out")[,,k])
	data6 = data6[,dim(data6)[2]:1]
	if (k == 3) data6 = data6/25.4
	if (k == 6) data6 = data6*100
	
	if (normalized) {
		
		start.date = year*10000 + 101
		end.date = year*10000 + 1231
		
		# NA nested grid:
		filename = "PM_geos5/data_conc_NA/dim.RData"
		date.vec = load.RData(filename, "date.vec")
		filename = paste("PM_geos5/data_conc_NA/conc_", spec.GC, ".RData", sep="")
		conc1 = apply(load.RData(filename, spec.GC)[,,which(date.vec >= start.date & date.vec <= end.date)], c(1,2), mean, na.rm=TRUE)
		
		# 2x2.5:
		filename = "PM_geos5/data_conc_2x2.5_201011/dim.RData"
		date.vec = load.RData(filename, "date.vec")
		filename = paste("PM_geos5/data_conc_2x2.5_201011/conc_", spec.GC, ".RData", sep="")
		conc2 = apply(load.RData(filename, spec.GC)[,,which(date.vec >= start.date & date.vec <= end.date)], c(1,2), mean, na.rm=TRUE)
		
		# 4x5:
		filename = "PM_geos5/data_conc_4x5_201011/dim.RData"
		date.vec = load.RData(filename, "date.vec")
		filename = paste("PM_geos5/data_conc_4x5_201011/conc_", spec.GC, ".RData", sep="")
		conc3 = apply(load.RData(filename, spec.GC)[,,which(date.vec >= start.date & date.vec <= end.date)], c(1,2), mean, na.rm=TRUE)
		
		# EPA-AQS / GEOS5:
		filename = "EPA_data/data_conc_AQS/AQS_GRID_2x2.5_2004-2008.RData"
		date.vec = load.RData(filename, "date.vec")
		if (species == "SOA" | species == "POA") specname = "OC"
		dataname = paste(specname, ".invdist", sep="")
		conc4 = apply(load.RData(filename, dataname)[,,which(date.vec >= start.date & date.vec <= end.date)], c(1,2), mean, na.rm=TRUE)
		
		# EPA-AQS / NCEP:
		if (species == "PM25") {
			filename = "EPA_data/data_conc_AQS/PM25_GRID_2.5x2.5_1998-2008a.RData"
			dataname = "PM25.day.invdist"
			}
		else {
			filename = "EPA_data/data_conc_AQS/SPEC_GRID_2.5x2.5_2000-2008.RData"
			if (species == "SOA" | species == "POA") specname = "OC"
			dataname = paste("PM25.", specname, ".invdist", sep="")
			}
		date.vec = load.RData(filename, "date.vec")
		conc5 = apply(load.RData(filename, dataname)[,dim(data5)[2]:1,which(date.vec >= start.date & date.vec <= end.date)], c(1,2), mean, na.rm=TRUE)
		
		# IMPROVE / NCEP:
		filename = "EPA_data/data_conc_IMPROVE/IMPROVE_GRID_2.5x2.5_1998-2006.RData"
		date.vec = load.RData(filename, "date.vec")
		if (species == "SO4") {
			dataname = "S.invdist"
			conc6 = apply((load.RData(filename, dataname)[,dim(data6)[2]:1,which(date.vec >= start.date & date.vec <= end.date)])*3, c(1,2), mean, na.rm=TRUE)
			}
		else {
			if (species == "SOA" | species == "POA") specname = "OC"
			dataname = paste(specname, ".invdist", sep="")
			conc6 = apply((load.RData(filename, dataname)[,dim(data6)[2]:1,which(date.vec >= start.date & date.vec <= end.date)]), c(1,2), mean, na.rm=TRUE)
			}
		
		data1 = data1/conc1*100
		data2 = data2/conc2*100
		data3 = data3/conc3*100
		data4 = data4/conc4*100
		data5 = data5/conc5*100
		data6 = data6/conc6*100
		
		}
	
	else { }
	
	# Decide zmax:
	
	if (is.null(zmax)) zmax = quantile(c(abs(as.vector(data1)), abs(as.vector(data2)), abs(as.vector(data3)), abs(as.vector(data4)), abs(as.vector(data5)), abs(as.vector(data6))), probs=0.975, na.rm=TRUE)
	
	data1[which(abs(data1) > zmax)] = sign(data1[which(abs(data1) > zmax)])*zmax
	data2[which(abs(data2) > zmax)] = sign(data2[which(abs(data2) > zmax)])*zmax
	data3[which(abs(data3) > zmax)] = sign(data3[which(abs(data3) > zmax)])*zmax
	data4[which(abs(data4) > zmax)] = sign(data4[which(abs(data4) > zmax)])*zmax
	data5[which(abs(data5) > zmax)] = sign(data5[which(abs(data5) > zmax)])*zmax
	data6[which(abs(data6) > zmax)] = sign(data6[which(abs(data6) > zmax)])*zmax
	
	# Plotting:
	
	if (normalized) {
		quartz(width=8, height=6.5, title=paste(species, "-", met, " Regression Coef (Normalized to Mean Conc in ", as.character(year), ")", sep=""))
		}
	
	else {
		quartz(width=8, height=6.5, title=paste(species, "-", met, " Regression Coef", sep=""))
		}
		
	par(mfrow=c(3,2))

	mai = c(0.1, 0.1, 0.2, 0.1)
	tcl = -0.1
	ps = 12
	legend.mar = 4.8

	filename = "PM_geos5/data_conc_NA/dim.RData"
	lon.plot = load.RData(filename, "lon.US")
	lat.plot = load.RData(filename, "lat.US")
	plot.field(data1, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar)
	add.text = "GEOS-Chem / GEOS5 - Nested NA"
	mtext(add.text, side=3)
	
	filename = "PM_geos5/data_conc_2x2.5_201011/dim.RData"
	lon.plot = load.RData(filename, "lon.US")
	lat.plot = load.RData(filename, "lat.US")
	plot.field(data2, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar)
	add.text = "GEOS-Chem / GEOS5 - 2x2.5"
	mtext(add.text, side=3)
	
	filename = "PM_geos5/data_conc_4x5_201011/dim.RData"
	lon.plot = load.RData(filename, "lon.US")
	lat.plot = load.RData(filename, "lat.US")
	plot.field(data3, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar)
	add.text = "GEOS-Chem / GEOS5 - 4x5"
	mtext(add.text, side=3)
	
	filename = "PM_geos5/data_conc_2x2.5_201011/dim.RData"
	lon.plot = load.RData(filename, "lon.US")
	lat.plot = load.RData(filename, "lat.US")
	plot.field(data4, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar)
	add.text = "EPA-AQS / GEOS5"
	mtext(add.text, side=3)
	
	filename = "EPA_data/data_conc_AQS/PM25_GRID_2.5x2.5_1998-2008a.RData"
	lon.plot = load.RData(filename, "lon.USmap")
	lat.plot = rev(load.RData(filename, "lat.US"))
	plot.field(data5, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar)
	add.text = "EPA-AQS / NCEP"
	mtext(add.text, side=3)
	
	filename = "EPA_data/data_conc_AQS/PM25_GRID_2.5x2.5_1998-2008a.RData"
	lon.plot = load.RData(filename, "lon.USmap")
	lat.plot = rev(load.RData(filename, "lat.US"))
	plot.field(data6, lon.plot, lat.plot, "sign", zlim=c(-zmax, zmax), mai=mai, tcl=tcl, ps=ps, legend.mar=legend.mar)
	add.text = "IMPROVE / NCEP"
	mtext(add.text, side=3)
	
	}



##################################################################


