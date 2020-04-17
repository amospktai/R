# This file contains useful functions for analyzing data obtained from the United States Environmental Protection Agency.
# Last update: Jan 2014, Amos P. K. Tai (pkamostai@gmail.com)



# Reading EPA monitoring site info from text file:
# This function reads the original monitoring site info text file, combines state code, county code and site ID into one site code, and reads the corresponding latitude and longitude.
# "file" has to have columns in this order: state code, county code, site ID, latitude, longitude.
get.site.info = function(filename, sep='\t', header=TRUE, na.strings=0) { 
    site.info = read.table(file=filename, sep=sep, header=header, na.strings=na.strings, colClasses=c(rep("character",3), rep("numeric",2)))
    site.info = cbind(paste(site.info[,1], site.info[,2], site.info[,3], sep=""), site.info[,4:5])
    site.info = na.omit(site.info)
    site.info[,1] = as.character(site.info[,1])
    colnames(site.info) = c("site", "lat", "lon")
    return(site.info)
}



# Reading generic site info from text or cvs file:
# This function reads the original monitoring site info text/cvs file with at least the site ID and the corresponding latitude and longitude.
# 'id.col', 'lat.col', 'lon.col' designate the column of the original data that contain info of site ID, latitude and longitude respectively.
get.site.info2 = function(filename, id.col, lat.col, lon.col, type='csv', sep=',', header=TRUE, na.strings=0) {
	if (type == 'csv') site.info = read.csv(file=filename, sep=sep, header=header, na.strings=na.strings) else site.info = read.table(file=filename, sep=sep, header=header, na.strings=na.strings, colClasses=c(rep("character",3), rep("numeric",2)))
	site.info = site.info[,c(id.col, lat.col, lon.col)]
	site.info = na.omit(site.info)
    site.info[,1] = as.character(site.info[,1])
    colnames(site.info) = c("site", "lat", "lon")
    return(site.info)
}



# Monthly statistics for data in date format "YYYYMMDD" (commonly used in EPA data):
# This function calculates monthly statistics for time series data of one year.
# "X" has two columns: 1. date in YYYYMMDD; 2. sample value.
# "fun" can be any basic sample statistics function e.g. "mean", "median", "max", etc.
monthstats3 <- function(X, year, fun) {
    monthly <- rep(0, times=12)
    ind1 <- which(X[,1] > (year*10000+100) & X[,1] < (year*10000+132))
    monthly[1] <- fun(X[ind1,2], na.rm=TRUE)  
    ind2 <- which(X[,1] > (year*10000+200) & X[,1] < (year*10000+232))
    monthly[2] <- fun(X[ind2,2], na.rm=TRUE)
    ind3 <- which(X[,1] > (year*10000+300) & X[,1] < (year*10000+332))
    monthly[3] <- fun(X[ind3,2], na.rm=TRUE)  
    ind4 <- which(X[,1] > (year*10000+400) & X[,1] < (year*10000+432))
    monthly[4] <- fun(X[ind4,2], na.rm=TRUE)
    ind5 <- which(X[,1] > (year*10000+500) & X[,1] < (year*10000+532))
    monthly[5] <- fun(X[ind5,2], na.rm=TRUE)
    ind6 <- which(X[,1] > (year*10000+600) & X[,1] < (year*10000+632))
    monthly[6] <- fun(X[ind6,2], na.rm=TRUE)
    ind7 <- which(X[,1] > (year*10000+700) & X[,1] < (year*10000+732))
    monthly[7] <- fun(X[ind7,2], na.rm=TRUE)
    ind8 <- which(X[,1] > (year*10000+800) & X[,1] < (year*10000+832))
    monthly[8] <- fun(X[ind8,2], na.rm=TRUE)
    ind9 <- which(X[,1] > (year*10000+900) & X[,1] < (year*10000+932))
    monthly[9] <- fun(X[ind9,2], na.rm=TRUE)
    ind10 <- which(X[,1] > (year*10000+1000) & X[,1] < (year*10000+1032))
    monthly[10] <- fun(X[ind10,2], na.rm=TRUE)
    ind11 <- which(X[,1] > (year*10000+1100) & X[,1] < (year*10000+1132))
    monthly[11] <- fun(X[ind11,2], na.rm=TRUE)
    ind12 <- which(X[,1] > (year*10000+1200) & X[,1] < (year*10000+1232))
    monthly[12] <- fun(X[ind12,2], na.rm=TRUE)
    return(monthly)
    }



# Monthly statistics for EPA AIRS data in date format "YYYYMMDD":
# This function extracts monthly statistics for one year of data directly from EPA AIRS data file ("file") in RD format, as downloaded from "http://www.epa.gov/ttn/airs/airsaqs/detaildata/downloadaqsdata.htm", and utilizes another monthly statistics function "monthstats3".
monthstats4 <- function(file, year, fun) {
    data <- read.table(file=file, sep="|", na.strings="", colClasses=c(rep("character",10), "numeric", "character", "numeric", rep("character",13), rep("numeric",2)))
    data <- cbind(paste(data[,3],data[,4],data[,5],sep=""), data[,c(6,9,11,13)])
    colnames(data) <- c("site", "parameter", "unit", "date", "sample.value")
    data <- na.omit(data)
    site.unique <- as.character(unique(data$site))
    monthly.all.site <- matrix(0, nrow=(length(site.unique)*12), ncol=3)
    monthly.all.site <- data.frame(monthly.all.site)
    for (i in 1:length(site.unique)) {
        ind <- which(data$site == site.unique[i])
        monthly.each.site <- monthstats3(data[ind, 4:5], year, fun)
        monthly.all.site[(i*12-11):(i*12),3] <- monthly.each.site   # col3 = sample value
        }
    monthly.all.site[,2] <- rep((year*10000+seq(1, 12)*100+1), times=length(site.unique))  # col2 = month
    monthly.all.site[,1] <- rep(site.unique, each=12) # col1 = site ID
    colnames(monthly.all.site) <- c("site", "month", "sample.value")
    return(monthly.all.site)
    }



# Reorganizing the monthly statistics for EPA AIRS data:
# This function reorganizes output "X" from "monthstats4" into another format with rows for distinct monitoring sites and columns for the twelve months.
monthstats5 <- function(X, year) {  # X is the output of monthstats4()
    site.unique <- unique(X[,1])
    M <- matrix(0, nrow=length(site.unique), ncol=12)
    for (i in 1:12) {
        M[,i] <- X[which(X[,2] == (year*10000+i*100+1)),3]
        }
    colnames(M) <- year*10000+seq(1,12)*100+1
    rownames(M) <- site.unique
    return(M)
    }



# Reading EPA AIRS data:
# This function reads EPA AIRS data file ("file") in RD format, as downloaded from "http://www.epa.gov/ttn/airs/airsaqs/detaildata/downloadaqsdata.htm".
# This function has the option (default) of removing all the rows where there are no sample values. The other fields can have "NA" as their entries.
read.EPA.RD = function(filename, na.rm=TRUE) {
	data = read.table(file=filename, sep="|", na.strings="", colClasses=c(rep("character",10), "numeric", "character", "numeric", rep("character",13), rep("numeric",2)))
    data = cbind(paste(data[,3], data[,4], data[,5],sep=""), data[,c(6:13,15,27,28)])
    colnames(data) = c("site", "parameter", "POC", "sample.duration", "unit", "method", "date", "time", "sample.value", "sample.freq", "detect.lim", "uncertainty")
    if (na.rm) data = data[which(!is.na(data$sample.value)),] else data = data
    return(data)
}      



# Reading EPA AIRS blank data:
# This function reads EPA AIRS blank data file ("file") in RB format, as downloaded from "http://www.epa.gov/ttn/airs/airsaqs/detaildata/downloadaqsdata.htm".
# This function removes all the rows where there are no sample values. The other fields can have "NA" as their entries.
read.EPA.RB = function(file) {
    data = read.table(file=file, sep="|", na.strings="", colClasses=c(rep("character",11), "numeric", "character", "numeric", rep("character",11), rep("numeric",2)))
    data = cbind(paste(data[,3],data[,4],data[,5],sep=""), data[,c(6,7,8,9,10,11,12,14,26,27)])
    colnames(data) <- c("site", "parameter", "POC", "sample.duration", "unit", "method", "blank.type", "date", "blank.value", "detect.lim", "uncertainty")
    data = data[which(is.na(data$blank.value) == FALSE),]
    return(data)
}      


	
# Reorganizing the daily statistics for EPA AIRS data:
# This function reorganizes the output ("data") from "read.EPA.RD" into another format with rows for distinct monitoring sites and columns for individual dates. When there are more than one sample for a given date, the mean is calculated.
read.EPA.RD2 = function(data) {
    site.unique = as.character(unique(data$site))
    date.unique = as.numeric(unique(data$date))
    date.unique = sort(date.unique)
    EPA.data = matrix(0, nrow=length(site.unique), ncol=length(date.unique))
    EPA.data = data.frame(EPA.data)
    for (i in 1:length(site.unique)) {
	ind.site = which(data$site == site.unique[i])
	data.subset = data[ind.site,]
	for (j in 1:length(date.unique)) {
		ind.date = which(data.subset$date == date.unique[j])
         	if (length(ind.date) == 1)
         		EPA.data[i,j] = data.subset$sample.value[ind.date]
         	else if (length(ind.date) > 1)
         	    EPA.data[i,j] = mean(data.subset$sample.value[ind.date], na.rm=TRUE)
         	else
         	   	EPA.data[i,j] = NaN
         	}
     	}
    rownames(EPA.data) = site.unique
    colnames(EPA.data) = date.unique
    return(EPA.data)
}


    
# Finding the latitude and longitude for a given monitoring site:
# This function finds the lat/lon values for a given monitoring by searching from a matrix of site info with columnes of: 1. site code; 2. latitude; 3. longitude. "site.info" can be output from function "get.site.info".
match.site.lat.lon = function(site.unique, site.info) {
    site.unique = as.character(site.unique)
    site.lat.lon = matrix(0, nrow=length(site.unique), ncol=2)
    for (i in 1:length(site.unique)) {
	site.ind = which(site.info[,1] == site.unique[i])
	if (length(site.ind) == 1) {
        site.lat.lon[i,1] = site.info[site.ind,2]
        site.lat.lon[i,2] = site.info[site.ind,3]
		}
	else
		site.lat.lon[i,] = NaN
        }
    colnames(site.lat.lon) = c("lat", "lon")
    return(site.lat.lon)
}


