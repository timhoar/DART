#======================================================================
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$
#
# Collection of functions used by DART
#
#======================================================================


#----------------------------------------------------------------------
# ncdf.make2 creates a netcdf file of some sort

ncdf.make2 <- function(dir = out.dir, name = "TEST", date = -1, ensemble, station ){

mv     <- -888888.
n.ens  <- length(ensemble)
n.stat <- length(station)
	
### define dimensions for scores

### counts from roc table
dim.roc <- dim.def.ncdf("Cont.columns", 
                        "Counts of Observed = TRUE, observed = FALSE", 1:2)

dim.rank <- dim.def.ncdf("rank.counts", "Counts of observed Ranks", 1:(n.ens + 1))

### number of forecasts of a specific prob, num obs, ntotal?
dim.rel   <- dim.def.ncdf("rel.counts", "Counts Prob. Used, Observed", 1:2)
dim.time  <- dim.def.ncdf("time", "Lead time", val = c(0,12,24,36,48,60) )
dim.level <- dim.def.ncdf("level", "hPa", 
             val = c(1000, 925,  850,  700,  500,  400,  300,  250,  200,  150,  100) )

dim.site  <- dim.def.ncdf("site", "id", val = 1:n.stat) ### hard coded - danger

###*** define date vector - reuse from original file?  unlimited?
dim.date    <- dim.def.ncdf("date",     "YYYYMMDDHH", date,     unlim = FALSE)
dim.ens     <- dim.def.ncdf("ensemble", "ensembles",  ensemble, unlim = FALSE)
dim.station <- dim.def.ncdf("station",  "station",    station,  unlim = FALSE)

###  thresholds for ROC plot
dim.thres   <- dim.def.ncdf("thresholds", "probability", 1:(n.ens + 1 ))

###**** one file per variable?  dimension for variable?  
windROC     <- var.def.ncdf("wind.roc", units = "int", 
               dim = list( dim.time, dim.date, dim.level, dim.thres, dim.roc),  
               missval = mv, 
               longname = "counts, based on exceeding a given threshold.")

windRANK    <- var.def.ncdf("wind.rank", units = "int", 
               dim = list(dim.time, dim.date, dim.level, dim.rank),
               missval = mv,
               longname = "Array with counts for rank histograms")

#windthres  <- var.def.ncdf("wind.thres", units = "m/s", 
#              dim = list(dim.level, dim.site),  
#              missval = mv, 
#              longname = "Threshold used to convert probs to binary values.")

wind.me     <- var.def.ncdf("wind.me", units = "m/s", 
               dim = list(dim.time, dim.date, dim.level),
               missval = mv,
               longname = "error (forecast - obs) for mean ensemble forecast") 

new <- create.ncdf(paste(out.dir, name, ".nc",sep = ""), 
                   list(windROC, windRANK, wind.me))
att.put.ncdf(new, 0, "source", name)
close.ncdf(new)
}

#----------------------------------------------------------------------
# calculate and remove bias 
#----------------------------------------------------------------------

bias.rm <- function(TT, bias.dims = c("ensemble","level","station")){

   ## Hacker et al: Tellus DOI:10.1111/j.1600-0870.2010.00497.x:
   ## 
   ##    "In a multiphysics or perturbed-parameter ensemble, each member
   ##    can be differently biased. Defining bias to be the experiment-mean 
   ##    error as a function of observing station, pressure level (or surface) 
   ##    and forecast lead time, we remove the bias of each forecast before
   ##    computing scores or spreads. That is, we use corrected individual-
   ##    member forecasts at a given lead time (f'_k = f_k - mean(o-f)),
   ##    where k indexes obervations or forecasts at a particular horizontal 
   ##    location, level and ensemble member, and the average is over 
   ##    available instances of verifying observations."
   ## 
   ## So, give me names of the dimensions to average over, (in any order)
   ## and I'll determine how to average the data array. 

   if(length(TT$dimstrings) != 5){ 
      stop("Bias removal is expecting an array with 5 dimensions")
   }

   ## Match the strings to the dimension indices - needed by 'apply'

   bias.id <- array(0,dim=length(TT$bias.dims))

   for (i in 1:length(TT$bias.dims)){
      bias.id[i] <- match(TT$bias.dims[i],TT$dimstrings)
      if (is.na(bias.id[i])){
         stop(sprintf("ERROR: cannot match bias string %s",TT$bias.dims[i]))
      }
   }

   ## seems more natural to have the ids in an ascending order, must
   ## permute the strings to match the ascending order.

   bob       <- sort(bias.id, index.return=T)
   bias.dims <- TT$bias.dims[bob$ix]
   bias.id   <- bias.id[  bob$ix]

   ## [obs - forecast] : positive bias => forecast too low.

   omf  <- TT$observed - TT$forecast
   bias <- apply(omf,   bias.id, mean, na.rm = TRUE)

   ## actually remove the bias. positive bias => forecast too low,
   ## must add the bias to the forecast.

   temp <- sweep(TT$forecast, bias.id, bias, FUN = "+" )

   ## The remaining quantities are not used AFAIK
   # mse <- apply(omf^2, bias.id, mean, na.rm = TRUE)
   # f   <- function(x){sum(is.finite(x), na.rm = TRUE)}
   # n   <- apply(omf, bias.id, f)

   return(list(data = temp, bias = bias))

}

#----------------------------------------------------------------------
# make output file - version 1
#----------------------------------------------------------------------
ncdf.make <- function(dir = out.dir, name = "TEST" ){
	
### define dimensions for scores

### counts from roc table
dim.roc   <- dim.def.ncdf("Cont.columns", 
                          "Counts of Observed = TRUE, observed = FALSE", 1:2)
dim.rank  <- dim.def.ncdf("rank.counts", "Counts of observed Ranks", 1:11)
dim.rel   <- dim.def.ncdf("rel.counts", "Counts Prob. Used, Observed", 1:2)
dim.time  <- dim.def.ncdf("time", "Lead time", val = c(0,12,24,36,48,60) )
dim.level <- dim.def.ncdf("level", "hPa", 
                   val = c(1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100) )
dim.site  <- dim.def.ncdf("site", "id", val = 1:31) ### hard coded - danger
date      <- dim.def.ncdf("date", "YYYYMMDDHH", -1, unlim = TRUE)

###  thresholds for ROC plot
threshold <- dim.def.ncdf("thresholds", "probability", (1:11 - 1)/10)

mv        <- -888888
windROC   <- var.def.ncdf("wind.roc", units = "int", 
                          dim = list(dim.roc, threshold, dim.level, dim.time, date),
                          missval = mv,
                          longname = "Array with counts, based on exceeding a given threshold.  ")

windRANK  <- var.def.ncdf("wind.rank", units = "int",
                          dim = list(dim.rank, dim.level, dim.time, date),
                          missval = mv, 
                          longname = "Array with counts for rank histograms")

windthres <- var.def.ncdf("wind.thres", units = "m/s",
                          dim = list(dim.level, dim.site),
                          missval = mv,
                          longname = "Threshold used to convert probs to binary values.")

wind.me   <- var.def.ncdf("wind.me", units = "m/s",
                          dim = list(dim.level, dim.site),
                          missval = mv,
                          longname = "error (forecast - obs) for mean ensemble forecast") 

new <- create.ncdf(paste(out.dir, name, ".nc",sep = ""),
                   list(windROC, windRANK, windthres, wind.me) )
att.put.ncdf( new, 0, "source", name)

close.ncdf(new)
}

#----------------------------------------------------------------------
# shape-preserving netCDF variable retrieval
# If there are singleton dimensions, they are preserved. Actually,
# they are reinstated. As such, it is pretty inefficient in that two 
# copies of the entire variable are required at the same time.
#
# It does return the dimension names, so you can check the storage
# order of the variable to make sure its what you expect.
# i.e. the 3rd dimension is 'level' or 'station' or 'longitude' ...
#----------------------------------------------------------------------
get.var.ncdf.whole <- function(nc, varname = NA, verbose = FALSE){

   if (class(nc) != "ncdf") {stop("nc not the result of open.ncdf()!")}

   # find the variable id that pertains to our variable name
   # and harvest all the dimension names and sizes

   varid <- -1 
   for (j in 1:nc$nvars){
       if ( nc$var[[j]]$name == varname ) {varid <- j}
   }
   if (varid < 1) {stop("varname not found")}

   naturalsize <- nc$var[[varid]]$size
   ndims       <- nc$var[[varid]]$ndims
   dimids      <- nc$var[[varid]]$dimids
   dims        <- nc$var[[varid]]$dim
   dimstring   <- character(ndims) 
   dimlens     <- integer(ndims)

   for (j in 1:ndims){
      dimstring[j] <- dims[[j]]$name 
      dimlens[  j] <- dims[[j]]$len 
   }

   vardata <- get.var.ncdf(nc,varname)
   alldata <- array(vardata,dim=naturalsize)

   return( list(file=nc$file, varname=varname, data=alldata, 
               dimstrings=dimstring, dimlengths=dimlens) )
}

#----------------------------------------------------------------------
# Does a variable exist in the netCDF file?
# Check without dying ...
#----------------------------------------------------------------------
exist.var.ncdf <- function(nc, varname = NA, verbose=FALSE){

   varid <- -1 

   if (class(nc) != "ncdf") {stop("nc not the result of open.ncdf()!")}

   # R actually keeps the coordinate variables separate from the 'variables'

   for (j in 1:nc$ndims){
       if ( nc$dim[[j]]$name == varname ) {varid <- j}
       if (verbose) print(sprintf('found a [%s] coordinate variable',nc$dim[[j]]$name))
   }

   # Check the variable names

   for (j in 1:nc$nvars){
       if ( nc$var[[j]]$name == varname ) {varid <- j}
       if (verbose) print(sprintf('found a [%s] variable',nc$var[[j]]$name))
   }

   if ((varid < 0) & verbose){
       warning(sprintf('Did not find a [%s] variable in %s',varname,nc$file))
   }

   return( varid )

}

#----------------------------------------------------------------------
# Search the dimstrings for desired dimnames and return the order. 
#----------------------------------------------------------------------
match.id.strings <- function(dimstrings=NULL, 
    mydimnames=c("time","date","ensemble","level","station","copy")){

   # Find dimension IDs - dimension names are different in 'old-style'
   # netCDF files than they are in 'modern' files.

   if (       is.finite(match('analysisT', dimstrings))) {
      id.analysis =     match('analysisT', dimstrings)
   } else if (is.finite(match('date',      dimstrings))) {
      id.analysis =     match('date',      dimstrings)
   } else {
      stop("Cannot find the analysis cycle time(s) in the netCDF file.")
   }

   if (       is.finite(match('forecast_lead', dimstrings))) {
      id.lead =         match('forecast_lead', dimstrings)
   } else if (is.finite(match('time',          dimstrings))) {
      id.lead =         match('time',          dimstrings)
   } else {
      stop("Cannot find the forecast lead times in the netCDF file.")
   }

   id.ens      = match('ensemble', dimstrings)
   id.level    = match('level',    dimstrings)
   id.station  = match('station',  dimstrings)
   id.copy     = match('copy',     dimstrings)

   return(list(lead  = id.lead,  ens     = id.ens,     copy     = id.copy, 
               level = id.level, station = id.station, analysis = id.analysis))

}
