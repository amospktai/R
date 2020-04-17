# This script takes CLM surface data and convert it into GC land cover input files. (Anthony Y.H. Wong, 2017)
# For details of using the 

# A function that flips cesm output map into our conventional GC mapping system
cesm_flip_output = function(map){
  temp <- map
  lon_length <- dim(map)[1]
  d <- length(dim(map))
  if (d == 1){
    map[1:(lon_length/2)] <- temp[(lon_length/2+1):lon_length] 
    map[(lon_length/2+1)] <- temp[1:(lon_length/2)]
  }
  if (d == 2){
    map[1:(lon_length/2),] <- temp[(lon_length/2+1):lon_length,] 
    map[(lon_length/2+1):lon_length,] <- temp[1:(lon_length/2),]
  }
  if (d == 3){
    map[1:(lon_length/2),,] <- temp[(lon_length/2+1):lon_length,,] 
    map[(lon_length/2+1):lon_length,,] <- temp[1:(lon_length/2),,]
  }
  if (d == 4){
    map[1:(lon_length/2),,,] <- temp[(lon_length/2+1):lon_length,,,] 
    map[(lon_length/2+1):lon_length,,,] <- temp[1:(lon_length/2),,,]
  }
  if (d == 5){
    map[1:(lon_length/2),,,,] <- temp[(lon_length/2+1):lon_length,,,,] 
    map[(lon_length/2+1):lon_length,,,,] <- temp[1:(lon_length/2),,,,]
  }
  return(map)
}


# Input File path. 
surfdata <- '../data_copy/surfdata_1.9x2.5_simyr2000_c130423.nc'

# Output file paths.
clm_landtype_file <- '../Future_Input/clm_landtype_soilcol_test.nc'
soilnox_file <- '../Future_Input/clm_soilnox_map_r85.nc'
aef_file <- '../Future_Input/aef_map_r85.nc'
pft_file <- '../Future_Input/pft_map_r85.nc'
lai_file <- '../Future_Input/lai_r45_new.nc'

#First create the deposition map
lake <- cesm_flip_output(load.nc(surfdata,'PCT_LAKE'))
wetland <- cesm_flip_output(load.nc(surfdata,'PCT_WETLAND'))
glacier <- cesm_flip_output(load.nc(surfdata,'PCT_GLACIER'))
urban <- apply(cesm_flip_output(load.nc(surfdata,'PCT_URBAN')),c(1,2),sum)
pft <- load.nc(surfdata,'PCT_PFT')

#Take out a particular segment of year
pft <- cesm_flip_output(pft)
landfrac <- cesm_flip_output(load.nc(surfdata,'LANDFRAC_PFT'))
area <- load.nc(surfdata,'AREA')
temp <- load.nc(surfdata,'LONGXY')[,1]
lat <- load.nc(surfdata,'LATIXY')[1,]
n <- length(temp)
lon <- 1:n
lon[1:(n/2)] <- temp[(n/2 +1):n] - 360
lon[(n/2 +1):n] <- temp[1:(n/2)]

clm_landtype <- array(0,c(dim(lake),22))
#lon <- seq(-180,177.5,2.5)
#lat <- seq(-90,90,length.out = 96)


clm_landtype[,,1] <- lake[,]/100
clm_landtype[,,19] <- glacier[,]/100
clm_landtype[,,20] <- urban[,]/100
clm_landtype[,,21] <- wetland[,]/100
clm_landtype[,,2:18] <- pft[,,]/100

#Handle landfrac
for (i in 1:22){
  clm_landtype[,,i] <- clm_landtype[,,i] * landfrac
}
clm_landtype[,,1] <- clm_landtype[,,1] + (1-landfrac)
area=area*1000^2

#Output 

x <- ncdim_def(name='lon', units='degrees_east', vals=lon, longname='longitude')
y <- ncdim_def(name='lat', units='degrees_north', vals=lat, longname='latitude')
lt <- ncdim_def(name='CLMLT',units='unitless',vals=1:22,longname = 'CLM land Types')
t <- ncdim_def(name='time', units='hours since 2005-01-01 00:00:00', vals=0, calendar='standard', longname='time')

var_olson <- ncvar_def(name='OLSON', units='unitless', dim=list(x, y, lt, t), 
                       longname='Native CLM gridded land type', prec='float', verbose=TRUE)
var_area = ncvar_def(name='DXYP', units='m2', dim=list(x, y, t), 
                     longname='Grid box surface area', prec='float', verbose=TRUE)

nc = nc_create(filename=clm_landtype_file, vars=list(var_olson,var_area,var_soilcol), force_v4 = TRUE, verbose=TRUE)
ncvar_put(nc, varid=var_olson, vals=clm_landtype, verbose=TRUE)
ncvar_put(nc, varid = var_area, vals = area, verbose = T)
ncvar_put(nc, varid = var_soilcol, vals = soil_color, verbose = T)

ncatt_put(nc, varid='lon', attname='axis', attval='X')
ncatt_put(nc, varid='lat', attname='axis', attval='Y')
ncatt_put(nc, varid='time', attname='axis', attval='T')
ncatt_put(nc, varid = 'CLMLT', attname = 'axis', attval = 'CLMLT')
ncatt_put(nc, varid=0, attname='Title', attval='DryDep Input Land Map from present day CLM surface data set')
ncatt_put(nc, varid=0, attname='Contact', attval='Anthony Y.H. Wong (ayhwong@cuhk.edu.hk)')
ncatt_put(nc, varid=0, attname='Conventions', attval='COARDS')
ncatt_put(nc, varid=0, attname='History', attval='15/12/2015: Preprocess with R')
nc_close(nc)

#Then create the soilnox map
# List of landtype:
# 1  Water
# 2  Collapsed
# 3  Ice
# 4  Cold Barren
# 5  Collapsed
# 6  Hot barren
# 7  collapsed
# 8  Open Shrubland ABC
# 9  Open Shrubland DE
# 10 Grassland/Savanna DE
# 11 Collapsed
# 12 Savannah/Grassland/Woody Savannah ABC
# 13 Collapsed
# 14 Collapsed
# 15 Mixed Forest
# 16 Evergreen broadleaf forest CDE
# 17 Deciduous broadleaf forest CDE
# 18 Deciduous needleleaf forest
# 19 Evergreen needleleaf forest
# 20 Deciduous broadleaf forest
# 21 Evergreen broadleaf forest
# 22 Crops
# 23 Urban
# 24 Collapsed

# Prepare 6 basic landtype
basic_landtype <- array(0,c(length(lon),length(lat),6))
basic_landtype[,,1] <- clm_landtype[,,1]                      # Water 
basic_landtype[,,2] <- clm_landtype[,,19]                     # Ice
basic_landtype[,,3] <- clm_landtype[,,20]                     # Urban
basic_landtype[,,4] <- clm_landtype[,,21]                     # Wetland
basic_landtype[,,5] <- clm_landtype[,,2]                      # Barren
basic_landtype[,,6] <- apply(clm_landtype[,,3:18],c(1,2),sum) # vegetation

# Prepare 4 basic vegetation type
vegtype <- array(0,c(length(lon),length(lat),4)) 
vegtype[,,1] <- apply(clm_landtype[,,3:10],c(1,2),sum)   # Tree
vegtype[,,2] <- apply(clm_landtype[,,11:13],c(1,2),sum)  # Shrub
vegtype[,,3] <- apply(clm_landtype[,,14:16],c(1,2),sum)  # Grass
vegtype[,,4] <- clm_landtype[,,17]                       # Crop

# Initialize Variable
soilnox_landtype <- array(0,c(length(lon),length(lat)))

#Main Loop
for (i in 1:length(lon)){
  for (j in 1:length(lat)){
    if (which.max(basic_landtype[i,j,]) == 1) {soilnox_landtype[i,j] <- 1}
    if (which.max(basic_landtype[i,j,]) == 2) {soilnox_landtype[i,j] <- 3}
    if (which.max(basic_landtype[i,j,]) == 3) {soilnox_landtype[i,j] <- 23}
    if (which.max(basic_landtype[i,j,]) == 4) {soilnox_landtype[i,j] <- 1}
    if (which.max(basic_landtype[i,j,]) == 5) {
      if (abs(lat[j]) > 60) {soilnox_landtype[i,j] <- 4} #distinguish cold and hot barren
      else {soilnox_landtype[i,j] <- 6}
    }
    # For vegetatation-dominated grid:
    if (which.max(basic_landtype[i,j,]) == 6){
      #For forest:
      if (which.max(vegtype[i,j,]) == 1){
        temp <- clm_landtype[i,j,3:10] / vegtype[i,j,1]
        tree_type <- c(temp[1]+temp[2],temp[3],temp[4],temp[5],temp[6],temp[7]+temp[8])
        # Assign mixed forest when no tree type is dominant
        if (max(tree_type) <= 0.7) {soilnox_landtype[i,j] <-15} 
        # If dominant tree type exist, assign accordingly
        if (tree_type[1] > 0.7) {soilnox_landtype[i,j] <- 19}
        if (tree_type[2] > 0.7) {soilnox_landtype[i,j] <- 18}
        if (tree_type[3] > 0.7) {soilnox_landtype[i,j] <- 21}
        if (tree_type[4] > 0.7) {soilnox_landtype[i,j] <- 16}
        if (tree_type[5] > 0.7) {soilnox_landtype[i,j] <- 20}
        if (tree_type[6] > 0.7) {soilnox_landtype[i,j] <- 17}
      }
      #For shrub
      if (which.max(vegtype[i,j,]) == 2){
        temp <- clm_landtype[i,j,11:13]
        if ((temp[1]+temp[2]) > temp[3]) {soilnox_landtype[i,j] <- 8}
        else {soilnox_landtype[i,j] <- 9}
      }
      #For grass
      if (which.max(vegtype[i,j,]) == 3){
        temp <- clm_landtype[i,j,14:16]
        if ((temp[2]+temp[3]) > temp[1]) {soilnox_landtype[i,j] <- 12}
        else {soilnox_landtype[i,j] <- 10}
      }
      #For crop
      if (which.max(vegtype[i,j,]) == 4) {soilnox_landtype[i,j] <- 22}
    }
  }
}

#Output to nc file 

x <- ncdim_def(name='lon', units='degrees_east', vals=lon, longname='longitude')
y <- ncdim_def(name='lat', units='degrees_north', vals=lat, longname='latitude')
t <- ncdim_def(name='time', units='hours since 2005-01-01 00:00:00', vals=0, calendar='standard', longname='time')

var_soilnox <- ncvar_def(name='LANDFRAC_K01', units='unitless', dim=list(x, y, t), 
                         longname='Fraction of land for each MODIS-KOPPEN landtype', 
                         prec='float', verbose=TRUE)
var_soilnox <- list(var_soilnox)

for (i in 2:24){
  temp <- ncvar_def(name=paste0('LANDFRAC_K',sprintf('%02d',i)), units='unitless', dim=list(x, y, t), 
                    longname='Fraction of land for each MODIS-KOPPEN landtype', 
                    prec='float', verbose=TRUE)
  var_soilnox <- append(var_soilnox,list(temp))
}

nc = nc_create(filename=soilnox_file, vars=var_soilnox, force_v4 = TRUE, verbose=TRUE)
for (i in 1:24){
  ncvar_put(nc, varid=var_soilnox[[i]], vals=(soilnox_landtype == i), verbose=TRUE)
}
ncatt_put(nc, varid='lon', attname='axis', attval='X')
ncatt_put(nc, varid='lat', attname='axis', attval='Y')
ncatt_put(nc, varid='time', attname='axis', attval='T')
ncatt_put(nc, varid=0, attname='Title', attval='SoilNOx Biome Map from CLM present day surface data')
ncatt_put(nc, varid=0, attname='Contact', attval='Anthony Y.H. Wong (ayhwong@cuhk.edu.hk)')
ncatt_put(nc, varid=0, attname='Conventions', attval='COARDS')
ncatt_put(nc, varid=0, attname='History', attval='31/8/2016: Made with R')
nc_close(nc)

#Then do the AEF and MEGAN PFT map

#Scaling factor 
f=10^-9/3600;

#PFT-Specific Data AEF
aef_isop=c(600,3000,1,7000,10000,7000,10000,11000,2000,4000, 4000, 1600, 800, 200, 1)*f;
aef_bpin=c(300 ,300 ,200 ,120 ,130 ,120 ,130 ,130 ,100, 150 ,100 ,1.5 ,1.5 ,1.5, 1.5)*f;
aef_care=c(160,160, 80 ,40 ,30, 40 ,30 ,30, 30 ,100, 30, 0.3 ,0.3 ,0.3 ,0.3)*f;
aef_limo=c(100,100, 130, 80, 80, 80, 80, 80, 60 ,100, 60, 0.7, 0.7, 0.7, 0.7)*f;
aef_ocim=c(70 ,70, 60 ,150 ,120 ,150, 120, 120, 90, 150, 90 ,2 ,2 ,2 ,2)*f;
aef_sabi=c(70,70, 40, 80, 50 ,80 ,50 ,50 ,50, 70 ,50 ,0.7, 0.7 ,0.7 ,0.7)*f;
aef_mbox=c(700, 60, 0.01, 0.01, 0.01, 0.01, 0.01, 2, 0.01, 0.01, 0.01, 0.01, 0.01 ,0.01, 0.01)*f

#Initialize variable
size=dim(clm_landtype);
aef_map_isop=array(0,dim = c(size[1],size[2]));
aef_map_bpin=array(0,dim = c(size[1],size[2]));
aef_map_care=array(0,dim = c(size[1],size[2]));
aef_map_limo=array(0,dim = c(size[1],size[2]));
aef_map_ocim=array(0,dim = c(size[1],size[2]));
aef_map_sabi=array(0,dim = c(size[1],size[2]));
aef_map_mbox=array(0,dim = c(size[1],size[2]))

#Start the multiplication
for (i in 1:15){
  aef_map_isop = aef_map_isop + clm_landtype [ , ,(i+2)] * aef_isop [i];
  aef_map_bpin = aef_map_bpin + clm_landtype [ , ,(i+2)] * aef_bpin [i];
  aef_map_care = aef_map_care + clm_landtype [ , ,(i+2)] * aef_care [i];
  aef_map_limo = aef_map_limo + clm_landtype [ , ,(i+2)] * aef_limo [i];
  aef_map_ocim = aef_map_ocim + clm_landtype [ , ,(i+2)] * aef_ocim [i];
  aef_map_sabi = aef_map_sabi + clm_landtype [ , ,(i+2)] * aef_sabi [i];
  aef_map_mbox = aef_map_mbox + clm_landtype [ , ,(i+2)] * aef_mbox [i];
}

#Write the AEF File 

x <- ncdim_def(name='lon', units='degrees_east', vals=lon, longname='longitude')
y <- ncdim_def(name='lat', units='degrees_north', vals=lat, longname='latitude')
t <- ncdim_def(name='time', units='hours since 2005-01-01 00:00:00', vals=0, calendar='standard', longname='time')

var_isop <- ncvar_def(name='AEF_ISOPRENE', units='kgC/m2/s', dim=list(x, y, t), 
                      longname='Emission Factor', prec='float', verbose=TRUE)
var_bpin <- ncvar_def(name='AEF_BETA_PINENE', units='kgC/m2/s', dim=list(x, y, t), 
                      longname='Emission Factor', prec='float', verbose=TRUE)
var_limo <- ncvar_def(name='AEF_LIMONENE', units='kgC/m2/s', dim=list(x, y, t), 
                      longname='Emission Factor', prec='float', verbose=TRUE)
var_care <- ncvar_def(name='AEF_CARENE', units='kgC/m2/s', dim=list(x, y, t), 
                      longname='Emission Factor', prec='float', verbose=TRUE)
var_ocim <- ncvar_def(name='AEF_OCIMENE', units='kgC/m2/s', dim=list(x, y, t), 
                      longname='Emission Factor', prec='float', verbose=TRUE)
var_sabi <- ncvar_def(name='AEF_SABINENE', units='kgC/m2/s', dim=list(x, y, t), 
                      longname='Emission Factor', prec='float', verbose=TRUE)
var_mbox <- ncvar_def(name='AEF_MBO', units='kgC/m2/s', dim=list(x, y, t), 
                      longname='Emission Factor', prec='float', verbose=TRUE)


nc = nc_create(filename=aef_file, 
               vars=list(var_isop,var_sabi,var_ocim,var_limo,var_care,var_bpin,var_mbox), 
               force_v4 = TRUE, verbose=TRUE)
ncvar_put(nc, varid = var_isop, vals = aef_map_isop, verbose=TRUE)
ncvar_put(nc, varid = var_bpin, vals = aef_map_bpin, verbose=TRUE)
ncvar_put(nc, varid = var_limo, vals = aef_map_limo, verbose=TRUE)
ncvar_put(nc, varid = var_care, vals = aef_map_care, verbose=TRUE)
ncvar_put(nc, varid = var_ocim, vals = aef_map_ocim, verbose=TRUE)
ncvar_put(nc, varid = var_sabi, vals = aef_map_sabi, verbose=TRUE)
ncvar_put(nc, varid = var_mbox, vals = aef_map_mbox, verbose=TRUE)

ncatt_put(nc, varid='lon', attname='axis', attval='X')
ncatt_put(nc, varid='lat', attname='axis', attval='Y')
ncatt_put(nc, varid='time', attname='axis', attval='T')
ncatt_put(nc, varid=0, attname='Title', attval='AEF Data set')
ncatt_put(nc, varid=0, attname='Contact', attval='Anthony Y.H. Wong (ayhwong@cuhk.edu.hk)')
ncatt_put(nc, varid=0, attname='Conventions', attval='COARDS')
ncatt_put(nc, varid=0, attname='History', attval='5/9/2016: Preprocess with R')
nc_close(nc)

x <- ncdim_def(name='lon', units='degrees_east', vals=lon, longname='longitude')
y <- ncdim_def(name='lat', units='degrees_north', vals=lat, longname='latitude')
t <- ncdim_def(name='time', units='hours since 2005-01-01 00:00:00', vals=0, calendar='standard', longname='time')

var1 <- ncvar_def(name='PFT_BARE', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var2 <- ncvar_def(name='PFT_NDLF_EVGN_TMPT_TREE', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var3 <- ncvar_def(name='PFT_NDLF_EVGN_BORL_TREE', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var4 <- ncvar_def(name='PFT_NDLF_DECD_BORL_TREE', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var5 <- ncvar_def(name='PFT_BDLF_EVGN_TROP_TREE', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var6 <- ncvar_def(name='PFT_BDLF_EVGN_TMPT_TREE', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var7 <- ncvar_def(name='PFT_BDLF_DECD_TROP_TREE', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var8 <- ncvar_def(name='PFT_BDLF_DECD_TMPT_TREE', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var9 <- ncvar_def(name='PFT_BDLF_DECD_BORL_TREE', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var10 <- ncvar_def(name='PFT_BDLF_EVGN_SHRB', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var11 <- ncvar_def(name='PFT_BDLF_DECD_TMPT_SHRB', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var12 <- ncvar_def(name='PFT_BDLF_DECD_BORL_SHRB', units='unitless', dim=list(x, y, t),
                  longname='PFT Coverage', prec='float', verbose=TRUE)
var13 <- ncvar_def(name='PFT_C3_ARCT_GRSS', units='unitless', dim=list(x, y, t),
                   longname='PFT Coverage', prec='float', verbose=TRUE)
var14 <- ncvar_def(name='PFT_C3_NARC_GRSS', units='unitless', dim=list(x, y, t),
                   longname='PFT Coverage', prec='float', verbose=TRUE)
var15 <- ncvar_def(name='PFT_C4_GRSS', units='unitless', dim=list(x, y, t),
                   longname='PFT Coverage', prec='float', verbose=TRUE)
var16 <- ncvar_def(name='PFT_CROP', units='unitless', dim=list(x, y, t),
                   longname='PFT Coverage', prec='float', verbose=TRUE)

nc = nc_create(filename=pft_file, 
               vars=list(var1,var2,var3,var4,var5,var6,var7,var8,
                         var9,var10,var11,var12,var13,var14,var15,var16), 
               force_v4 = TRUE, verbose=TRUE)


ncvar_put(nc,varid = "PFT_BARE",vals = clm_landtype[ , ,2]);
ncvar_put(nc,varid = "PFT_NDLF_EVGN_TMPT_TREE",vals = clm_landtype[ , ,3]);
ncvar_put(nc,varid = "PFT_NDLF_EVGN_BORL_TREE",vals = clm_landtype[ , ,4]);
ncvar_put(nc,varid = "PFT_NDLF_DECD_BORL_TREE",vals = clm_landtype[ , ,5]);
ncvar_put(nc,varid = "PFT_BDLF_EVGN_TROP_TREE",vals = clm_landtype[ , ,6]);
ncvar_put(nc,varid = "PFT_BDLF_EVGN_TMPT_TREE",vals = clm_landtype[ , ,7]);
ncvar_put(nc,varid = "PFT_BDLF_DECD_TROP_TREE",vals = clm_landtype[ , ,8]);
ncvar_put(nc,varid = "PFT_BDLF_DECD_TMPT_TREE",vals = clm_landtype[ , ,9]);
ncvar_put(nc,varid = "PFT_BDLF_DECD_BORL_TREE",vals = clm_landtype[ , ,10]);
ncvar_put(nc,varid = "PFT_BDLF_EVGN_SHRB",vals = clm_landtype[ , ,11]);
ncvar_put(nc,varid = "PFT_BDLF_DECD_TMPT_SHRB",vals = clm_landtype[ , ,12]);
ncvar_put(nc,varid = "PFT_BDLF_DECD_BORL_SHRB",vals = clm_landtype[ , ,13]);
ncvar_put(nc,varid = "PFT_C3_ARCT_GRSS",vals = clm_landtype[ , ,14]);
ncvar_put(nc,varid = "PFT_C3_NARC_GRSS",vals = clm_landtype[ , ,15]);
ncvar_put(nc,varid = "PFT_C4_GRSS",vals = clm_landtype[ , ,16]);
ncvar_put(nc,varid = "PFT_CROP",vals = clm_landtype[ , ,17]);

ncatt_put(nc, varid='lon', attname='axis', attval='X')
ncatt_put(nc, varid='lat', attname='axis', attval='Y')
ncatt_put(nc, varid='time', attname='axis', attval='T')
ncatt_put(nc, varid=0, attname='Title', attval='AEF Data set')
ncatt_put(nc, varid=0, attname='Contact', attval='Anthony Y.H. Wong (ayhwong@cuhk.edu.hk)')
ncatt_put(nc, varid=0, attname='Conventions', attval='COARDS')
ncatt_put(nc, varid=0, attname='History', attval='5/9/2016: Preprocess with R')
nc_close(nc)

nc_close(nc);

#Finally, LAI
pftlai <- load.nc(surfdata,'MONTHLY_LAI')

pftlai <- cesm_flip_output(pftlai)
pftlai_out <- array(0,c(144,96,22,12))
pftlai_out[,,2:18,] <- pftlai

#Output 

x <- ncdim_def(name='lon', units='degrees_east', vals=lon, longname='longitude')
y <- ncdim_def(name='lat', units='degrees_north', vals=lat, longname='latitude')
lt <- ncdim_def(name='CLMLT',units='unitless',vals=1:22,longname = 'CLM land Types')
t <- ncdim_def(name='time', units='Months since 2005-01-01 00:00:00',vals=1:12, calendar='standard', longname='time')

var_modis <- ncvar_def(name='MODIS', units='unitless', dim=list(x, y, lt, t), 
                       longname='PFT-specific LAI', prec='float', verbose=TRUE)

nc = nc_create(filename=lai_file, vars=var_modis, force_v4 = TRUE, verbose=TRUE)
ncvar_put(nc, varid=var_modis, vals=pftlai_out, verbose=TRUE)

ncatt_put(nc, varid='lon', attname='axis', attval='X')
ncatt_put(nc, varid='lat', attname='axis', attval='Y')
ncatt_put(nc, varid='time', attname='axis', attval='T')
ncatt_put(nc, varid = 'CLMLT', attname = 'axis', attval = 'CLMLT')
ncatt_put(nc, varid=0, attname='Title', attval='PFT-Specific LAI')
ncatt_put(nc, varid=0, attname='Contact', attval='Anthony Y.H. Wong (ayhwong@cuhk.edu.hk)')
ncatt_put(nc, varid=0, attname='Conventions', attval='COARDS')
ncatt_put(nc, varid=0, attname='History', attval='15/12/2015: Preprocess with R')
nc_close(nc)
