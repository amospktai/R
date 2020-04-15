# Load GEOS-Chem data:
load.ctm = function(PATH='./', CAT='IJ-AVG-$', TRC=1, YEAR, MONTH, DAY=1, HOUR=0, LON=c(-180, 180), LAT=c(-90, 90), LEV=NULL, FUN=NULL, apply.over=3) {
   YYYYMMDDhh = NULL
   for (y in 1:length(YEAR)) {
      YYYY = as.character(YEAR[y])
      for (m in 1:length(MONTH)) {
         MM = as.character(MONTH[m])
         if (MONTH[m] < 10) MM = paste('0', MM, sep='')
         for (d in 1:length(DAY)) {
            DD = as.character(DAY[d])
            if (DAY[d] < 10) DD = paste('0', DD, sep='')
            for (h in 1:length(HOUR)) {
               hh = as.character(HOUR[h])
               if (HOUR[h] < 10) hh = paste('0', hh, sep='')
               YYYYMMDDhh = c(YYYYMMDDhh, paste(YYYY, MM, DD, hh, sep=''))
            }
         }
      }
   }
   TTT = as.character(TRC)
   if (TRC < 10) TTT = paste('00', TTT, sep='') else if (TRC >= 10 & TRC < 100) TTT = paste('0', TTT, sep='') else TTT = TTT
   # Get full dimensions:
   filename = paste(PATH, 'ctm_', chartr(old='$', new='S', CAT), '_N', TTT, '_', YYYYMMDDhh[1], '.nc', sep='')
   lon = load.nc(filename, varid='lon')
   lat = load.nc(filename, varid='lat')
   lev = load.nc(filename, varid='lev')
   # Extract entire data:
   X = array(NaN, dim=c(length(lon), length(lat), length(lev), length(YYYYMMDDhh)))
   for (t in 1:length(YYYYMMDDhh)) {
      filename = paste(PATH, 'ctm_', chartr(old='$', new='S', CAT), '_N', TTT, '_', YYYYMMDDhh[t], '.nc', sep='')
      X[,,,t] = load.nc(filename, varid=paste(CAT, '_N', TTT, sep=''))
   }
   # Extract subset if any:
   if (is.null(LEV)) LEV = lev
   ind_lon = which(lon >= LON[1] & lon <= LON[2])
   ind_lat = which(lat >= LAT[1] & lat <= LAT[2])
   X = X[ind_lon,ind_lat,LEV,]
   LON = lon[ind_lon]
   LAT = lat[ind_lat]
   if (!is.null(FUN)) X = apply(X, MARGIN=(1:length(dim(X)))[-apply.over], FUN=FUN, na.rm=TRUE)
   out = list(var=X, lon=LON, lat=LAT, lev=LEV, time=YYYYMMDDhh)
   return(out)
}
