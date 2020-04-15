##################################################################
##################################################################



# Map showing observed vs. simulated concentrations at native and 2x2.5 resolutions:

plot.conc.compare.mapx3 = function(species, year=2006) {
	
	# Plot observations, native-resolution and 2x2.5 concentrations.
	
	# Use Cooke OC emissions, non-local PBL mixing, no SOA dependence on inorganic aerosol mass.
	
	# Mean site-by-site concentrations for 2006 only:
	site.info = get.site.info(file="~/Dropbox/Research/EPA_data/EPA_site_latlon.txt", sep="\t", header=TRUE, na.strings=0)

	if (species == "PM25") spec.GC = "PMtotchm"
	if (species == "SO4") spec.GC = "sulfate"
	if (species == "NO3") spec.GC = "nitrate"
	if (species == "OC") spec.GC = "OC"
	if (species == "EC") spec.GC = "BC"
	
	start.date = year*10000 + 101
	end.date = year*10000 + 1231
	
	# Native resolution:
	filename = '~/Dropbox/Research/PM_geos5/data_conc_NA/dim.RData'
	lon.NA = load.RData(filename, "lon.US")
	lat.NA = load.RData(filename, "lat.US")
	date.NA = load.RData(filename, "date.vec")
	filename = paste('~/Dropbox/Research/PM_geos5/data_conc_NA/conc_', spec.GC, '.RData', sep='')
	data.NA = load.RData(filename, spec.GC)[,,which(date.NA >= start.date & date.NA <= end.date)]
	plot.NA = apply(data.NA, c(1,2), mean, na.rm=TRUE)
	
	if (species == 'PM25') {
		# Need to scale down OC:
		date.X = date.NA
		filename = '~/Dropbox/Research/PM_geos5/data_conc_NA/conc_OC.RData'
		data.X = load.RData(filename, 'OC')[,,which(date.X >= start.date & date.X <= end.date)]
		OC.adjust = apply(data.X, c(1,2), mean, na.rm=TRUE)
		plot.NA = plot.NA - OC.adjust*1.4*(1 - 0.9029)
	}
		
	if (species == 'OC') {
		# Need to scale down OC to remove dependence on inorganic aerosol:
		plot.NA = plot.NA*0.9029
	}
	
	# 2x2.5 resolution:
	filename = '~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/dim.RData'
	lon.2x2.5 = load.RData(filename, "lon.US")
	lat.2x2.5 = load.RData(filename, "lat.US")
	date.2x2.5 = load.RData(filename, "date.vec")
	filename = paste('~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/conc_', spec.GC, '.RData', sep='')
	data.2x2.5 = load.RData(filename, spec.GC)[,,which(date.2x2.5 >= start.date & date.2x2.5 <= end.date)]
	plot.2x2.5 = apply(data.2x2.5, c(1,2), mean, na.rm=TRUE)
	
	if (species == 'PM25') {
		# Need to subtract 201011 OC and add 201008 OC:
		filename = '~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/dim.RData'
		date.X = load.RData(filename, 'date.vec')
		filename = '~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/conc_OC.RData'
		data.X = load.RData(filename, 'OC')[,,which(date.X >= start.date & date.X <= end.date)]
		OC.subtract = apply(data.X, c(1,2), mean, na.rm=TRUE)
		plot.2x2.5 = plot.2x2.5 - OC.subtract*1.4
		filename = '~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201008/dim.RData'
		date.X = load.RData(filename, 'date.vec')
		filename = '~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201008/conc_OC.RData'
		data.X = load.RData(filename, 'OC')[,,which(date.X >= start.date & date.X <= end.date)]
		OC.add = apply(data.X, c(1,2), mean, na.rm=TRUE)
		# Scale up from full mixing to non-local mixing:
		plot.2x2.5 = plot.2x2.5 + (OC.add*1.051 - 0.052)*1.4
	}
		
	if (species == 'OC') {
		# Use 201008 OC instead and scale from full mixing to non-local mixing:
		filename = '~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201008/dim.RData'
		date.X = load.RData(filename, 'date.vec')
		filename = '~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201008/conc_OC.RData'
		data.X = load.RData(filename, 'OC')[,,which(date.X >= start.date & date.X <= end.date)]
		plot.2x2.5 = apply(data.X, c(1,2), mean, na.rm=TRUE)
		plot.2x2.5 = plot.2x2.5*1.051 - 0.052
	}
	
	# Observations:
	
	find.site.avg = function(data.in) {
		site = as.vector(unique(data.in$site))
		data.out = matrix(0, nrow=length(site), ncol=3)
		for (n in 1:length(site)) {
			data.out[n,1] = mean(data.in$sample.value[which(data.in$site == site[n])], na.rm=TRUE)
			data.out[n,2:3] = match.site.lat.lon(site[n], site.info)
		}
		return(data.out)
	}
		
	find.site.adj = function(site.avg, blank.avg) {
		site.avg.adj = site.avg
		for (i in 1:nrow(site.avg.adj)) {
			ind.colocated = which(blank.avg[,2] == site.avg[i,2] & blank.avg[,3] == site.avg[i,3])
			colocated = length(ind.colocated)
			if (colocated == 0)	site.avg.adj[i,2:3] = NaN
			else site.avg.adj[i,1] = site.avg[i,1] - blank.avg[ind.colocated,1]
			if (site.avg.adj[i,1] < 0) site.avg.adj[i,1] = NaN
		}
		return(site.avg.adj)
	}
	
	filename = "~/Dropbox/Research/EPA_data/data_conc_AQS/AQS_SITE_2004-2008.RData"
	data.obs = load.RData(filename, paste(species, ".data", sep=""))
	data.obs.sub = data.obs[which(data.obs$date >= start.date & data.obs$date <= end.date),]
	date.obs.sub = data.obs$date[which(data.obs$date >= start.date & data.obs$date <= end.date)]
			
	site.avg = find.site.avg(data.obs.sub)
			
	if (species == "OC") {
			
	blank.data = load.RData("~/Dropbox/Research/EPA_data/data_conc_AQS/OC_BLANK.RData", "OC.blank.data")
	data.in = blank.data[which(blank.data$date >= start.date & blank.data$date <= end.date),]
	site = as.vector(unique(data.in$site))
	data.out = matrix(0, nrow=length(site), ncol=3)
	for (n in 1:length(site)) {
		data.out[n,1] = mean(data.in$blank.value[which(data.in$site == site[n])], na.rm=TRUE)
		data.out[n,2:3] = match.site.lat.lon(site[n], site.info)
	}
	blank.avg = data.out
	site.avg = find.site.adj(site.avg, blank.avg)
			
	}
			
	else { }
	
	zmax = max(c(as.vector(plot.NA), as.vector(plot.2x2.5), site.avg[,1]), na.rm=TRUE)
		
	mai = c(0.05, 0.15, 0.05, 0.15)
	mgp = c(0.8, 0.4, 0)
	tcl = -0.2
	ps = 16
	legend.mar = 2.5
	legend.width = 2.0
	
	plot.field(plot.NA, lon.NA, lat.NA, type="def", zlim=c(0, zmax), mai=mai, mgp=mgp, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
	
	plot.field(plot.2x2.5, lon.2x2.5, lat.2x2.5, type="def", zlim=c(0, zmax), mai=mai, mgp=mgp, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
	
	plot.site(site.avg[,1], site.avg[,2:3], xlim=c(head(lon.NA, 1), tail(lon.NA, 1)), ylim=c(head(lat.NA, 1), tail(lat.NA, 1)), type="def", zlim=c(0, zmax), pch=20, mai=mai, mgp=mgp, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
	#text(x=-122, y=26, labels=species, cex=2)
	
	OUTPUT = list(conc.NA=plot.NA, conc.2x2.5=plot.2x2.5, conc.obs=site.avg)
	return(OUTPUT)
		
}



##################################################################



# Scatter plots for observed vs. simulated concentrations or concentrations at different resolutions:

plot.conc.compare.scat = function(resol1, conc.1, resol2=NULL, conc.2=NULL, obs=FALSE, conc.obs=NULL, method='spavg', divide.EW=FALSE) {
	
	# conc.1 is a 3-d array with 3rd dim representing species 'SO4', 'NO3' and 'OC'.
	# conc.obs is a list of species 'SO4', 'NO3' and 'OC'.
	
	filename = '~/Dropbox/Research/PM_geos5/data_conc_NA/dim.RData'
	lon.NA = load.RData(filename, "lon.US")
	lat.NA = load.RData(filename, "lat.US")
	filename = '~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/dim.RData'
	lon.2x2.5 = load.RData(filename, "lon.US")
	lat.2x2.5 = load.RData(filename, "lat.US")
	filename = '~/Dropbox/Research/PM_geos5/data_conc_4x5_201011/dim.RData'
	lon.4x5 = load.RData(filename, "lon.US")
	lat.4x5 = load.RData(filename, "lat.US")

	if (!is.null(resol2)) {
		
		if (resol1 == 'NA' & resol2 == '2x2.5') {
			SO4.NA = conc.1[,,1]
			NO3.NA = conc.1[,,2]
			OC.NA = conc.1[,,3]
			SO4.2x2.5 = conc.2[,,1]
			NO3.2x2.5 = conc.2[,,2]
			OC.2x2.5 = conc.2[,,3]
			SO4.X = as.vector(sp.dissolve(SO4.NA[,length(lat.NA):1], lon.NA, rev(lat.NA), lon.2x2.5, rev(lat.2x2.5))[,length(lat.2x2.5):1])
			NO3.X = as.vector(sp.dissolve(NO3.NA[,length(lat.NA):1], lon.NA, rev(lat.NA), lon.2x2.5, rev(lat.2x2.5))[,length(lat.2x2.5):1])
			OC.X = as.vector(sp.dissolve(OC.NA[,length(lat.NA):1], lon.NA, rev(lat.NA), lon.2x2.5, rev(lat.2x2.5))[,length(lat.2x2.5):1])
			SO4.Y = as.vector(SO4.2x2.5)
			NO3.Y = as.vector(NO3.2x2.5)
			OC.Y = as.vector(OC.2x2.5)
		}
	
		if (resol1 == 'NA' & resol2 == '4x5') {
			SO4.NA = conc.1[,,1]
			NO3.NA = conc.1[,,2]
			OC.NA = conc.1[,,3]
			SO4.4x5 = conc.2[,,1]
			NO3.4x5 = conc.2[,,2]
			OC.4x5 = conc.2[,,3]
			SO4.X = as.vector(sp.dissolve(SO4.NA[,length(lat.NA):1], lon.NA, rev(lat.NA), lon.4x5, rev(lat.4x5))[,length(lat.4x5):1])
			NO3.X = as.vector(sp.dissolve(NO3.NA[,length(lat.NA):1], lon.NA, rev(lat.NA), lon.4x5, rev(lat.4x5))[,length(lat.4x5):1])
			OC.X = as.vector(sp.dissolve(OC.NA[,length(lat.NA):1], lon.NA, rev(lat.NA), lon.4x5, rev(lat.4x5))[,length(lat.4x5):1])
			SO4.Y = as.vector(SO4.4x5)
			NO3.Y = as.vector(NO3.4x5)
			OC.Y = as.vector(OC.4x5)
		}
	
		if (resol1 == '2x2.5' & resol2 == '4x5') {
			SO4.2x2.5 = conc.1[,,1]
			NO3.2x2.5 = conc.1[,,2]
			OC.2x2.5 = conc.1[,,3]
			SO4.4x5 = conc.2[,,1]
			NO3.4x5 = conc.2[,,2]
			OC.4x5 = conc.2[,,3]
			SO4.X = as.vector(sp.dissolve(SO4.2x2.5[,length(lat.2x2.5):1], lon.2x2.5, rev(lat.2x2.5), lon.4x5, rev(lat.4x5))[,length(lat.4x5):1])
			NO3.X = as.vector(sp.dissolve(NO3.2x2.5[,length(lat.2x2.5):1], lon.2x2.5, rev(lat.2x2.5), lon.4x5, rev(lat.4x5))[,length(lat.4x5):1])
			OC.X = as.vector(sp.dissolve(OC.2x2.5[,length(lat.2x2.5):1], lon.2x2.5, rev(lat.2x2.5), lon.4x5, rev(lat.4x5))[,length(lat.4x5):1])
			SO4.Y = as.vector(SO4.4x5)
			NO3.Y = as.vector(NO3.4x5)
			OC.Y = as.vector(OC.4x5)
		}

		reg.SO4 = lmodel2(SO4.Y ~ SO4.X)
		reg.NO3 = lmodel2(NO3.Y ~ NO3.X)
		reg.OC = lmodel2(OC.Y ~ OC.X)
		
		for (x in 1:3) {
			
			if (x == 1) {
				X = SO4.X; Y = SO4.Y; reg = reg.SO4
				X.name = expression(S*O[4]^'2-')
			}
			
			if (x == 2) {
				X = NO3.X; Y = NO3.Y; reg = reg.NO3
				X.name = expression(N*O[3]^'-')
			}
			
			if (x == 3) {
				X = OC.X; Y = OC.Y; reg = reg.OC
				X.name = 'OC'
			}

			r.sq = reg$rsquare
			a = reg$regression.results[2,2]
			b = reg$regression.results[2,3]
			#NMB = sum((Y - X), na.rm=TRUE)/sum(X, na.rm=TRUE)
			NMB = sum((X - Y), na.rm=TRUE)/sum(Y, na.rm=TRUE)
			xylim = c(0, max(c(X, Y)))
			plot(X, Y, xlim=xylim, ylim=xylim, xlab='', ylab='')
			abline(a=0, b=1)
			abline(a=a, b=b, lty=2)
			text(x=(xylim[1] + 0.082*diff(xylim)), y=(xylim[1] + 0.961*diff(xylim)), labels=expression(R^2*' ='), cex=1.6)
			text(x=(xylim[1] + 0.268*diff(xylim)), y=(xylim[1] + 0.950*diff(xylim)), labels=as.character(round(r.sq, digits=1.5)), cex=1.6)
			text(x=(xylim[1] + 0.220*diff(xylim)), y=(xylim[1] + 0.850*diff(xylim)), labels=paste('b = ', as.character(round(b, digits=3)), sep=''), cex=1.6)
			text(x=(xylim[1] + 0.220*diff(xylim)), y=(xylim[1] + 0.750*diff(xylim)), labels=paste('NMB = ', as.character(round(NMB, digits=3)), sep=''), cex=1.6)
			text(x=(xylim[1] + 0.85*diff(xylim)), y=(xylim[1] + 0.15*diff(xylim)), labels=X.name, cex=2.0)
			
		}

	}
	
	if (is.null(resol2) & obs) {
		
		SO4.obs = conc.obs[[1]]
		NO3.obs = conc.obs[[2]]
		OC.obs = conc.obs[[3]]
				
		if (resol1 == 'NA') {
			SO4.NA = conc.1[,,1]
			NO3.NA = conc.1[,,2]
			OC.NA = conc.1[,,3]
			SO4.res = SO4.NA
			NO3.res = NO3.NA
			OC.res = OC.NA
			lon.res = lon.NA
			lat.res = lat.NA
		}
	
		if (resol1 == '2x2.5') {
			SO4.2x2.5 = conc.1[,,1]
			NO3.2x2.5 = conc.1[,,2]
			OC.2x2.5 = conc.1[,,3]
			SO4.res = SO4.2x2.5
			NO3.res = NO3.2x2.5
			OC.res = OC.2x2.5
			lon.res = lon.2x2.5
			lat.res = lat.2x2.5
		}
	
		if (resol1 == '4x5') {
			SO4.4x5 = conc.1[,,1]
			NO3.4x5 = conc.1[,,2]
			OC.4x5 = conc.1[,,3]
			SO4.res = SO4.4x5
			NO3.res = NO3.4x5
			OC.res = OC.4x5
			lon.res = lon.4x5
			lat.res = lat.4x5
		}
		
		mai=c(0.2, 0.2, 0.05, 0.05)
		mgp=c(1.0, 0.4, 0)
		cex=1.3

		reg.SO4 = plot.mod.obs.US(SO4.res, SO4.obs, lon.res, lat.res, xlab='', ylab='', mai=mai, mgp=mgp, cex=cex, method=method, divide.EW=divide.EW)
		reg.NO3 = plot.mod.obs.US(NO3.res, NO3.obs, lon.res, lat.res, xlab='', ylab='', mai=mai, mgp=mgp, cex=cex, method=method, divide.EW=divide.EW)
		reg.OC = plot.mod.obs.US(OC.res, OC.obs, lon.res, lat.res, xlab='', ylab='', mai=mai, mgp=mgp, cex=cex, method=method, divide.EW=divide.EW)
		
	}
	
	LIST = list(SO4=reg.SO4, NO3=reg.NO3, OC=reg.OC)
	return(LIST)
	
}



##################################################################



# Map showing observed vs. simulated standard deviations in 2x2.5:

plot.conc.sd.compare.mapx2 = function(species, mai=c(0.05,0.15,0.05,0.15), mgp=c(0.8,0.4,0), tcl=-0.2, ps=16, legend.mar=2.5, legend.width=2.0, model.left=TRUE) {
	
	if (species == "PM25") spec.GC = "PMtotchm"
	if (species == "SO4") spec.GC = "sulfate"
	if (species == "NO3") spec.GC = "nitrate"
	if (species == "OC") spec.GC = "OC"
	if (species == "EC") spec.GC = "BC"
	
	filename = '~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/dim.RData'
	lon.US = load.RData(filename, "lon.US")
	lat.US = load.RData(filename, "lat.US")
	
	# Simulated:
	filename = paste('~/Dropbox/Research/PM_geos5/data_conc_2x2.5_201011/conc_', spec.GC, '.RData', sep='')
	conc.std = load.RData(filename, paste(spec.GC, '.std', sep=''))
	data1 = rm.nonUS.2x2.5(apply(conc.std, c(1,2), sd, na.rm=TRUE))
	
	# Observed:
	filename = '~/Dropbox/Research/EPA_data/data_conc_AQS/AQS_GRID_2x2.5_2004-2008.RData'
	conc.std = load.RData(filename, paste(species, '.std', sep=''))
	data2 = rm.nonUS.2x2.5(apply(conc.std, c(1,2), sd, na.rm=TRUE))
	
	zmax = max(c(as.vector(data1), as.vector(data2)), na.rm=TRUE)
	
	if (model.left) {
		
		plot.field(data1, lon.US, lat.US, type="def", zlim=c(0, zmax), mai=mai, mgp=mgp, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
	
		plot.field(data2, lon.US, lat.US, type="def", zlim=c(0, zmax), mai=mai, mgp=mgp, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
		
	} else {
		
		plot.field(data2, lon.US, lat.US, type="def", zlim=c(0, zmax), mai=mai, mgp=mgp, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
		
		plot.field(data1, lon.US, lat.US, type="def", zlim=c(0, zmax), mai=mai, mgp=mgp, tcl=tcl, ps=ps, legend.mar=legend.mar, legend.width=legend.width)
		
	}
	
	OUTPUT = list(mod=data1, obs=data2)
	return(OUTPUT)
	
}



##################################################################
##################################################################



# Plot function for comparison:

plot.conc.compare.season = function(species, obs.data='AQS', mod.data='2x2.5_201011', start.date=20060101, end.date=20061231, method="spavg") {
	
	# Mean site-by-site concentrations for 2006 only:
site.info = get.site.info(file="~/Dropbox/Research/EPA_data/EPA_site_latlon.txt", sep="\t", header=TRUE, na.strings=0)

	# Modeled concentrations:
	
	species = species
	if (species == "PM25") spec.GC = "PMtotchm"
	if (species == "SO4") spec.GC = "sulfate"
	if (species == "NO3") spec.GC = "nitrate"
	if (species == "OC") spec.GC = "OC"
	if (species == "EC") spec.GC = "BC"
	
	filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", mod.data, "/dim.RData", sep="")
	lon.mod = load.RData(filename, "lon.US")
	lat.mod = load.RData(filename, "lat.US")
	date.mod = load.RData(filename, "date.vec")
	date.mod.sub = date.mod[which(date.mod >= start.date & date.mod <= end.date)]
	
	filename = paste("~/Dropbox/Research/PM_geos5/data_conc_", mod.data, "/conc_", spec.GC, ".RData", sep="")
	data.mod = load.RData(filename, spec.GC)
	data.mod.sub = data.mod[,,which(date.mod >= start.date & date.mod <= end.date)]
	plot.mod.ann = apply(data.mod.sub, c(1,2), mean, na.rm=TRUE)
	plot.mod.DJF = find.season.mean(data.mod.sub, date.mod.sub, "DJF")
	plot.mod.MAM = find.season.mean(data.mod.sub, date.mod.sub, "MAM")
	plot.mod.JJA = find.season.mean(data.mod.sub, date.mod.sub, "JJA")
	plot.mod.SON = find.season.mean(data.mod.sub, date.mod.sub, "SON")
	
	
	# Observed concentrations:
	
	find.site.avg = function(data.in) {
		site = as.vector(unique(data.in$site))
		data.out = matrix(0, nrow=length(site), ncol=3)
		for (n in 1:length(site)) {
			data.out[n,1] = mean(data.in$sample.value[which(data.in$site == site[n])], na.rm=TRUE)
			data.out[n,2:3] = match.site.lat.lon(site[n], site.info)
			}
		return(data.out)
		}
		
	find.site.adj = function(site.avg, blank.avg) {
		site.avg.adj = site.avg
		for (i in 1:nrow(site.avg.adj)) {
			ind.colocated = which(blank.avg[,2] == site.avg[i,2] & blank.avg[,3] == site.avg[i,3])
			colocated = length(ind.colocated)
			if (colocated == 0)	site.avg.adj[i,2:3] = NaN
			else site.avg.adj[i,1] = site.avg[i,1] - blank.avg[ind.colocated,1]
			if (site.avg.adj[i,1] < 0) site.avg.adj[i,1] = NaN
			}
		return(site.avg.adj)
		}
	
	{
	
	if (obs.data == "AQS" | obs.data == "AQS+IMPROVE") {
		
		filename = "~/Dropbox/Research/EPA_data/data_conc_AQS/AQS_SITE_2004-2008.RData"
		data.obs = load.RData(filename, paste(species, ".data", sep=""))
		
		data.obs.sub = data.obs[which(data.obs$date >= start.date & data.obs$date <= end.date),]
		date.obs.sub = data.obs$date[which(data.obs$date >= start.date & data.obs$date <= end.date)]
		date.obs.sub.new = date.obs.sub - signif(date.obs.sub, digits=4)
		
		site.avg.ann = find.site.avg(data.obs.sub)
		site.avg.DJF = find.site.avg(data.obs.sub[which(date.obs.sub.new <= 229 | date.obs.sub.new >= 1201),])
		site.avg.MAM = find.site.avg(data.obs.sub[which(date.obs.sub.new >= 301 & date.obs.sub.new <= 531),])
		site.avg.JJA = find.site.avg(data.obs.sub[which(date.obs.sub.new >= 601 & date.obs.sub.new <= 831),])
		site.avg.SON = find.site.avg(data.obs.sub[which(date.obs.sub.new >= 901 & date.obs.sub.new <= 1131),])
		
		{
		
		if (species == "OC") {
			
			blank.data = load.RData("~/Dropbox/Research/EPA_data/data_conc_AQS/OC_BLANK.RData", "OC.blank.data")
			data.in = blank.data[which(blank.data$date >= start.date & blank.data$date <= end.date),]
			site = as.vector(unique(data.in$site))
			data.out = matrix(0, nrow=length(site), ncol=3)
			for (n in 1:length(site)) {
				data.out[n,1] = mean(data.in$blank.value[which(data.in$site == site[n])], na.rm=TRUE)
				data.out[n,2:3] = match.site.lat.lon(site[n], site.info)
				}
			blank.avg = data.out

			site.avg.ann = find.site.adj(site.avg.ann, blank.avg)
			site.avg.DJF = find.site.adj(site.avg.DJF, blank.avg)
			site.avg.MAM = find.site.adj(site.avg.MAM, blank.avg)
			site.avg.JJA = find.site.adj(site.avg.JJA, blank.avg)
			site.avg.SON = find.site.adj(site.avg.SON, blank.avg)
			
			}
			
		else { }
		
		}
		
		{
		
		if (obs.data == "AQS+IMPROVE") {
			
			filename = "~/Dropbox/Research/EPA_data/data_conc_IMPROVE/IMPROVE_SITE_1998-2006.RData"
			if (species == "SO4") {
				data.obs = load.RData(filename, "S.data")
				data.obs$sample.value = data.obs$sample.value*3
				}
			else {
				data.obs = load.RData(filename, paste(species, ".data", sep=""))
				}
		
			data.obs.sub = data.obs[which(data.obs$date >= start.date & data.obs$date <= end.date),]
			date.obs.sub = data.obs$date[which(data.obs$date >= start.date & data.obs$date <= end.date)]
			date.obs.sub.new = date.obs.sub - signif(date.obs.sub, digits=4)
		
			site.avg.ann.2 = find.site.avg(data.obs.sub)
			site.avg.DJF.2 = find.site.avg(data.obs.sub[which(date.obs.sub.new <= 229 | date.obs.sub.new >= 1201),])
			site.avg.MAM.2 = find.site.avg(data.obs.sub[which(date.obs.sub.new >= 301 & date.obs.sub.new <= 531),])
			site.avg.JJA.2 = find.site.avg(data.obs.sub[which(date.obs.sub.new >= 601 & date.obs.sub.new <= 831),])
			site.avg.SON.2 = find.site.avg(data.obs.sub[which(date.obs.sub.new >= 901 & date.obs.sub.new <= 1131),])
			
			site.avg.ann = rbind(site.avg.ann, site.avg.ann.2)
			site.avg.DJF = rbind(site.avg.DJF, site.avg.DJF.2)
			site.avg.MAM = rbind(site.avg.MAM, site.avg.MAM.2)
			site.avg.JJA = rbind(site.avg.JJA, site.avg.JJA.2)
			site.avg.SON = rbind(site.avg.SON, site.avg.SON.2)
			
			}
		
		else { }
		
		}
		
		}
		
	else if (obs.data == "IMPROVE") {
		
		filename = "~/Dropbox/Research/EPA_data/data_conc_IMPROVE/IMPROVE_SITE_1998-2006.RData"
		if (species == "SO4") {
			data.obs = load.RData(filename, "S.data")
			data.obs$sample.value = data.obs$sample.value*3
			}
		else {
			data.obs = load.RData(filename, paste(species, ".data", sep=""))
			}
		
		data.obs.sub = data.obs[which(data.obs$date >= start.date & data.obs$date <= end.date),]
		date.obs.sub = data.obs$date[which(data.obs$date >= start.date & data.obs$date <= end.date)]
		date.obs.sub.new = date.obs.sub - signif(date.obs.sub, digits=4)
		
		site.avg.ann = find.site.avg(data.obs.sub)
		site.avg.DJF = find.site.avg(data.obs.sub[which(date.obs.sub.new <= 229 | date.obs.sub.new >= 1201),])
		site.avg.MAM = find.site.avg(data.obs.sub[which(date.obs.sub.new >= 301 & date.obs.sub.new <= 531),])
		site.avg.JJA = find.site.avg(data.obs.sub[which(date.obs.sub.new >= 601 & date.obs.sub.new <= 831),])
		site.avg.SON = find.site.avg(data.obs.sub[which(date.obs.sub.new >= 901 & date.obs.sub.new <= 1131),])
		
		}
	
	else {
		print("There is no corresponding observed dataset.")
		stop()
		}
			
	}
		
	quartz(title=paste("Mod-Obs Comparison (GEOS-Chem ", mod.data, " vs. ", obs.data, ") for ", species, sep=""), width=8, height=7)
	par(mfrow=c(5, 3))
	
	zmax = max(c(as.vector(plot.mod.ann), as.vector(plot.mod.DJF), as.vector(plot.mod.MAM), as.vector(plot.mod.JJA), as.vector(plot.mod.SON), site.avg.ann[,1], site.avg.DJF[,1], site.avg.MAM[,1], site.avg.JJA[,1], site.avg.SON[,1]), na.rm=TRUE)
	
	for (t in 1:5) {
		
		if (t == 1) {
			plot.mod = plot.mod.ann
			site.avg = site.avg.ann
			lab = "Annual"
			}
		else if (t == 2) {
			plot.mod = plot.mod.DJF
			site.avg = site.avg.DJF
			lab = "DJF"
			}
		else if (t == 3) {
			plot.mod = plot.mod.MAM
			site.avg = site.avg.MAM
			lab = "MAM"
			}
		else if (t == 4) {
			plot.mod = plot.mod.JJA
			site.avg = site.avg.JJA
			lab = "JJA"
			}
		else {
			plot.mod = plot.mod.SON
			site.avg = site.avg.SON
			lab = "SON"
			}
		
		mai = c(0.2, 0.2, 0.1, 0.1)
		mgp = c(0.7, 0.2, 0)
		tcl = -0.2
		ps = 10
		legend.mar = 2
		
		plot.site(site.avg[,1], site.avg[,2:3], xlim=c(head(lon.mod, 1), tail(lon.mod, 1)), ylim=c(head(lat.mod, 1), tail(lat.mod, 1)), type="def", zlim=c(0, zmax), pch=20, mai=mai, mgp=mgp, tcl=tcl, ps=ps, legend.mar=legend.mar)
		
		plot.field(plot.mod, lon.mod, lat.mod, type="def", zlim=c(0, zmax), mai=mai, mgp=mgp, tcl=tcl, ps=ps, legend.mar=legend.mar)
	
		mai = c(0.2, 0.3, 0.1, 1.1)
		mgp = c(0.7, 0.15, 0)
		tcl = -0.1
		ps = 10
		
		plot.mod.obs.US(mod.data=plot.mod, site.data=site.avg, lon.mod=lon.mod, lat.mod=lat.mod, method=method, n=2, dlim=500, mai=mai, mgp=mgp, tcl=tcl, ps=ps)
		mtext(text=lab, side=4, line=1.5, cex=2)
		
		}
		
	}



##################################################################
##################################################################


