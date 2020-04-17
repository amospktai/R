################################################################################
### Examples of using "ncdf4" and "rhdf5" packages in R
################################################################################

################################################################################
### EXAMPLE 1: (ncdf4)
### Prepare GEOS-Chem NH3 base emission files without agricultural sources
################################################################################

# Objective: Extract MASAGE agricultural NH3 and EMEP total NH3 emission -> subtract MASAGE agricultural NH3 from EMEP total NH3 to get EMEP non-agricultural NH3 -> store it as nc file

setwd('/scratch/prjtgabi/gcgrid/data/ExtData/HEMCO')
load('/scratch/prjtgabi/R/functions.RData')
library(ncdf4); library(fields); library(maps)


# MASAGE agricultural NH3 emission:

MASAGE_var = c('COTTON', 'FRUITS_VEGETABLES', 'GROUNDNUT', 'MAIZE', 'OIL_PALM', 'OTHER_CEREALS', 'OTHER_CEREALS_WINTER', 'OTHER_CROPS', 'PULSES', 'RAPESEED', 'EARLY_RICE', 'LATE_RICE', 'SOYBEANS', 'SUGAR_CANE', 'TEMPERATE_ROOTS', 'TROPICAL_CEREALS', 'TROPICAL_ROOTS', 'WHEAT', 'WHEAT_WINTER', 'BEEF_CROP', 'BEEF_HOUSING', 'BEEF_PASTURE', 'BEEF_STORAGE', 'BUFFALOES_EMISSION', 'DAIRY_CROP', 'DAIRY_HOUSING', 'DAIRY_PASTURE', 'DAIRY_STORAGE', 'GOATS_EMISSION', 'PORK_CROP', 'PORK_HOUSING', 'PORK_PASTURE', 'PORK_STORAGE', 'POULTRY_CROP', 'POULTRY_HOUSING', 'POULTRY_PASTURE', 'POULTRY_STORAGE', 'SHEEP_EMISSION')

# Extract data from nc file:
nc = nc_open('MASAGE_NH3/v2015-02/masage_nh3_2006.2x25.nc'); nc
lon_2x25 = ncvar_get(nc, 'lon')
lat_2x25 = ncvar_get(nc, 'lat')
time_MASAGE = ncvar_get(nc, 'time')
mon_MASAGE = c("2006-01-01", "2006-02-01", "2006-03-01", "2006-04-01", "2006-05-01", "2006-06-01", "2006-07-01", "2006-08-01", "2006-09-01", "2006-10-01", "2006-11-01", "2006-12-01")
MASAGE_NH3 = array(NaN, dim=c(length(lon_2x25), length(lat_2x25), length(time_MASAGE), length(MASAGE_var)))
for (l in 1:length(MASAGE_var)) MASAGE_NH3[,,,l] = ncvar_get(nc, MASAGE_var[l])
nc_close(nc)

# Calculate total (over all subtypes) and annual mean NH3:
MASAGE_NH3_total = apply(MASAGE_NH3, 1:3, sum, na.rm=TRUE)
MASAGE_NH3_mean = apply(MASAGE_NH3_total, 1:2, mean, na.rm=TRUE)


# EMEP total NH3 emission:

# Extract data from nc file:
nc = nc_open('EMEP/v2015-03/EMEP.generic.1x1.nc')
lon_EMEP = ncvar_get(nc, 'lon')
lat_EMEP = ncvar_get(nc, 'lat')
time_EMEP = ncvar_get(nc, 'time')
EMEP_NH3 = ncvar_get(nc, 'NH3')
nc_close(nc)

# Find fraction of agricultural sources:
ind_2006 = 17
EMEP_NH3_2006_2x25 = sp.regrid(EMEP_NH3[,,ind_2006], lon.in=lon_EMEP, lat.in=lat_EMEP, lon.out=lon_2x25, lat.out=lat_2x25)
frac_agr = MASAGE_NH3_mean/EMEP_NH3_2006_2x25
frac_agr[which(frac_agr == 'Inf')] = 1
frac_agr[which(frac_agr > 1)] = 1
plot.field(frac_agr[61:93,61:81], lon_2x25[61:93], lat_2x25[61:81])
frac_agr_EMEP = sp.regrid(frac_agr, lon.in=lon_2x25, lat.in=lat_2x25, lon.out=lon_EMEP, lat.out=lat_EMEP)
frac_agr_EMEP[which(is.na(frac_agr_EMEP))] = 0
plot.field(frac_agr_EMEP, lon_EMEP, lat_EMEP)

# Create new variable for non-agricultural NH3:
EMEP_NH3_nagr = array(NaN, dim=dim(EMEP_NH3))
for (k in 1:dim(EMEP_NH3)[3]) EMEP_NH3_nagr[,,k] = EMEP_NH3[,,k]*(1 - frac_agr_EMEP)

# Create new nc file:
# Dimensions:
x = ncdim_def(name='lon', units='degrees_east', vals=lon_EMEP, longname='longitude')
y = ncdim_def(name='lat', units='degrees_north', vals=lat_EMEP, longname='latitude')
t = ncdim_def(name='time', units='hours since 1990-01-01 00:00:00', vals=time_EMEP, calendar='standard', longname='time')
# Variables:
var_NH3 = ncvar_def(name='NH3', units='kg/m2/s', dim=list(x, y, t), longname='NH3 tracer (non-agricultural)', prec='float', verbose=TRUE)
nc = nc_create(filename='EMEP/v2015-03/EMEP.NH3.nagr.1x1.nc', vars=var_NH3, force_v4=TRUE, verbose=TRUE)
# If multiple variables are to be put into the same file, "vars" should be a list of objects returned by "ncvar_def", e.g., vars=list(var_NH3, var_NH4, var_NO3).
# Putting in values:
ncvar_put(nc, varid=var_NH3, vals=EMEP_NH3_nagr, verbose=TRUE)
ncatt_put(nc, varid='lon', attname='axis', attval='X')
ncatt_put(nc, varid='lat', attname='axis', attval='Y')
ncatt_put(nc, varid='time', attname='axis', attval='T')
ncatt_put(nc, varid=0, attname='Title', attval='EMEP non-agricultural NH3')
ncatt_put(nc, varid=0, attname='Contact', attval='Amos P. K. Tai (amostai@cuhk.edu.hk)')
ncatt_put(nc, varid=0, attname='Conventions', attval='COARDS')
ncatt_put(nc, varid=0, attname='History', attval='May 2015: subtract MASAGE agricultural NH3 from EMEP total NH3')
nc_close(nc)



################################################################################
### EXAMPLE 2: (rhdf5)
### Extract and analyze satellite-derived PM2.5 data in China
################################################################################

# Objective: Extract satellite-derived annual mean PM2.5 from hdf5 files, regrid them into coarser resolution and save them as RData files.


setwd('/users/staff/b132099/projects/PM_clim_CN/')
load('/users/staff/b132099/R/functions.RData')
library(rhdf5); library(fields); library(maps)

filename1 = 'Unified_PM25_CH_1998_to_2012/Unified_PM25_CH_'
filename2 = '01_'
filename3 = '12-RH50-minc0_Median.h5'
year_name = as.character(seq(1998, 2012))

y = 1
filename = paste(filename1, year_name[y], filename2, year_name[y], filename3, sep='')
lon_CN = h5read(filename, 'LON')
lat_CN = h5read(filename, 'LAT')

PM25_unfilter = array(NaN, dim=c(length(lon_CN), length(lat_CN), length(year_name)))
for (y in 1:length(year_name)) {
	filename = paste(filename1, year_name[y], filename2, year_name[y], filename3, sep='')
	PM25_unfilter[,,y] = t(h5read(filename, 'PM25'))[,length(lat_CN):1]/10
}
PM25_CN = PM25_unfilter

# Remove negative values:
PM25_CN[which(PM25_CN <= 0)] = NaN

# Where are the outliers?
ind_outlier = which(PM25_CN > 300, arr.ind=TRUE)
PM25_outlier = data.frame(cbind(lon_CN[ind_outlier[,1]], lat_CN[ind_outlier[,2]], as.numeric(year_name[ind_outlier[,3]]), PM25_CN[ind_outlier]))
colnames(PM25_outlier) = c('lon', 'lat', 'year', 'PM25')
write.csv(PM25_outlier, file='PM25_outlier.csv', row.names=FALSE)

# Remove outliers (>99.99% quantile or >200ug/m3):
# PM25_CN[which(PM25_CN > quantile(PM25_CN, 0.9999, na.rm=TRUE))] = NaN
PM25_CN[which(PM25_CN > 200)] = NaN

# China region for NCEP/NCAR native grid:
# lon: 72.5 (ind = 30), 135.0 (ind = 55)
# lat: 17.5 (ind = 30), 55.0 (ind = 15)

lon = seq(72.5, 135.0, by=2.5)
lat = seq(17.5, 55.0, by=2.5)

# Regrid from 0.1x0.1 to 2.5x2.5:
PM25 = sp.dissolve.3D(PM25_CN, lon_CN, lat_CN, lon, lat)

y = 3
#quartz.par(mfrow=c(1,2), width=10, height=4)
plot.field(PM25_CN[,,y], lon_CN, lat_CN)
plot.field(PM25[,,y], lon, lat)

# Save data:
save(list=c('PM25', 'lon', 'lat', 'year_name'), file='PM25_CN_annual_2.5x2.5_1998-2012.RData')
save(list=c('PM25_CN', 'lon_CN', 'lat_CN', 'year_name'), file='PM25_CN_annual_0.1x0.1_1998-2012.RData')


######################################################################
######################################################################

