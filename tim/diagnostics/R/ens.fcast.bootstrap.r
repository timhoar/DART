ens.fcast.bootstrap <- function(filename=NULL, varname=NULL, quantiles=NULL, 
                        add.obs.error=TRUE, bias.dims=NULL, group.dims=NULL,
                        n.boot=100, prob.method='quantile', fixed.thres=NA, 
                        run.bs=TRUE, run.crps=TRUE, run.rank=TRUE, run.attr=TRUE, 
                        out.file=NULL) {
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
# DART:get.var.ncdf.whole() maintains singleton array dimensions
# and returns a list with the storage order. The required storage order is 
#
# in netcdf: (analysisT, stations, levels, copy, nmembers, forecast_lead)
# in      R:  forecast_lead, nmembers, copy, levels, stations, analysisT 
# 
# Since the analysisT dimension is the record/unlimited dimension, multiple
# forecast experiments may be concatenated into a single netCDF file.
# Previously, each forecast exeriment (each analysis cycle) was contained
# in a separate file and it was necessary to loop over them ...
#======================================================================

#----------------------------------------------------------------------
# READ INPUT
#----------------------------------------------------------------------

nc        <- open.ncdf(filename)
stations  <- get.var.ncdf(nc, "station")
level     <- get.var.ncdf(nc, "level")
copy      <- get.var.ncdf(nc, "copy")
ensemble  <- get.var.ncdf(nc, "ensemble")

# TJH FIXME units for these
# The obs_seq_verify netCDF files use 'forecast_lead',
# the old netCDF files use 'time' 
if (         exist.var.ncdf(nc,"forecast_lead") > 0){
   leads    <- get.var.ncdf(nc,"forecast_lead")
   timfriendly <- TRUE
} else if (  exist.var.ncdf(nc,"time") > 0){
   leads    <- get.var.ncdf(nc,"time")*3600   ## convert hours to seconds
   timfriendly <- FALSE
} else {
   stop('nothing good can come from this')
}

# TJH FIXME units for these
# The obs_seq_verify netCDF files use 'analysisT',
# the old netCDF files use 'date' 
if (          exist.var.ncdf(nc,"analysisT") > 0){
   analysisT <- get.var.ncdf(nc,"analysisT")
} else if (   exist.var.ncdf(nc,"date") > 0){
   analysisT <- get.var.ncdf(nc,"date")
} else {
   stop('nothing good can come from this either')
}

copy.meta <- sub(' +$','',get.var.ncdf(nc, "CopyMetaData")) # trim trailing whitespace

n.analysis <- length(analysisT)
n.station  <- length(stations)
n.level    <- length(level)
n.ens      <- length(ensemble)
n.lead     <- length(leads)
n.copy     <- length(copy)

# Determine the indices pertaining to the
# forecast, observation, and observation error variance.

if ( timfriendly ) {
   obindex = match('observation',copy.meta)
   oeindex = match('observation error variance',copy.meta)
   fcindex = match('forecast',copy.meta)
} else {
   obindex = match('observations (metar)',copy.meta)
   oeindex = match('observation error standard deviation',copy.meta)
   fcindex = match('observation - forecast',copy.meta)
}

if (is.finite(obindex)){ 
   print(sprintf('observation       copy is index %d',obindex))
} else {
   stop("Cannot find the observation copy index.")
}

if (is.finite(fcindex)){ 
   print(sprintf('forecast          copy is index %d',fcindex))
} else {
   stop("Cannot find the forecast copy index.")
}

if (is.finite(oeindex)){ 
   print(sprintf('observation error copy is index %d',oeindex))
} else {
   stop("Cannot find the observation error copy index.")
}
   
# If the netCDF file has a global attribute short name, use it to identify
# the experiment.

NMS <- att.get.ncdf(nc, 0, "SHORT_NAME")$value
print(NMS)

#----------------------------------------------------------------------
# ncdf.functions.r:get.var.ncdf.whole() maintains singleton array 
# dimensions and returns a list with the storage order.
# We read this and permute it to the following storage order:
#
# R: [fcst_lead, ensemble, copy, level, station, analysisT]
#
# TJH - I'd really like to remove the copy dimension and have :
#
# W$forecast [n.lead, n.ens, n.level, n.station, n.analysis]
# W$observed [n.lead, n.ens, n.level, n.station, n.analysis]
# W$oberrvar [n.lead, n.ens, n.level, n.station, n.analysis]
#
# and then in the bootstrap loop (depends on the set of analyses)
#
# boot$forecast [n.lead, n.ens, n.level, n.station, n.analysis]
# boot$observed [n.lead,     1, n.level, n.station, n.analysis]
# boot$oberrvar [n.lead,     1, n.level, n.station, n.analysis]
#----------------------------------------------------------------------

# v.units <- att.get.ncdf(nc, varname, "units"    )$value
long.name <- att.get.ncdf(nc, varname, "long_name")$value
missvalue <- att.get.ncdf(nc, varname, "missing_value")$value

W         <- get.var.ncdf.whole(nc, varname)

# Find dimension IDs - dimension names are different in 'old-style'
# netCDF files than they are in 'modern' files.

id <- match.id.strings(W$dimstrings)

if ( is.na(id$lead)   | is.na(id$ens)      | is.na(id$copy) |
     is.na(id$level)  | is.na(id$station)  | is.na(id$analysis)){
   stop("Cannot find the required dimension names in the netCDF file.")
}

W$dimstrings[id$lead]     = 'forecast_lead'
W$dimstrings[id$analysis] = 'analysisT'

#----------------------------------------------------------------------
# permute the dimensions to a known, preferred shape
# R: [fcst_lead, ensemble, copy, level, station, analysisT]
#----------------------------------------------------------------------

myorder <- c(id$lead, id$ens, id$copy, id$level, id$station, id$analysis)

W$data       <- aperm(W$data,perm=myorder)
W$dimstrings <- W$dimstrings[myorder]
W$dimlengths <- W$dimlengths[myorder]
id           <- match.id.strings(W$dimstrings) # recompute ids given new order

if ( ! timfriendly ) {

   # These netCDF files did not correctly use the missing_value
   # attribute - they incorrectly used the _FillValue ...
   # They also had the observation error standard deviation (not variance)
   # and they had (O-F), not the actual F ...

   W$data[W$data < (1+missvalue)] <- NA # replace missing values 

   ## stored observation error standard deviation
   ## need observation error variance

   W$data[,,oeindex,,,] = W$data[,,oeindex,,,,drop=F]^2

   copy.meta[oeindex] <- 'observation error variance'

   ## actually stored O-F in forecast copy slot ... reconstruct forecast
   ## observation - (observation - forecast) == forecast

   W$data[,,fcindex,,,] = W$data[,,obindex,,,,drop=F] - 
                          W$data[,,fcindex,,,,drop=F]

   copy.meta[fcindex] <- 'forecast'

}

## Extract forecast/observed/oberrvar and remove the 'copy' dimension.
## For the sake of the bias removal, just keep identical shapes
## for observed and oberrvar for now.
## forecast  [lead, ensemble, level, station, analysisT]
## observed  [lead, ensemble, level, station, analysisT]

dims <- W$dimlengths[-id$copy] 

org = list( forecast = array(W$data[,,fcindex,,,],dim=dims),
            observed = array(W$data[,,obindex,,,],dim=dims),
            oberrvar = array(W$data[,,oeindex,,,],dim=dims),
            dimstrings = W$dimstrings[-id$copy],
            file=W$file, varname=W$varname)

id <- match.id.strings(org$dimstrings) # recompute ids given new order

rm(W)

neworder <- c(id$lead, id$ens, id$level, id$station, id$analysis)

# TJH ... losing the whole one-level-at-a-time philosophy.
# Doing all levels, so I commented out the following block.
#
# TJH FIXME - P should be between 1 and n.level
# Use a single level and maintain singleton dimensions
#
# if ( W$dimstrings[id$level] == "level" ) {
#    print(sprintf('Using level %f (of %d possible)\n',P[1],n.level))
#    bootdata <- W$data[ , , ,P[1], ,drop=F]
# } else {
#    stop(sprintf("dimension %d is not 'level'",id$level))
# }

# TJH DEBUG
# scatterplot of all forecast leads for all stations - obs vs. (O-F)
# for the first ensemble member, first level, first analysis
plot(org$observed[,1,1,,1],org$forecast[,1,1,,1], main=long.name, 
       sub=filename, xlab=copy.meta[obindex], ylab=copy.meta[fcindex])

# This block left as a guideline to demonstrate how wind speed was
# created and processed originally. Unsupported/Untested in this code.
# different storage order ... 

if(varname == "WS"){
   V <- get.var.ncdf.whole(nc, "V")
   U <- get.var.ncdf.whole(nc, "U")

   if(length(dim(U)) != 6){ stop("wrong dimensions!") }

   U <- U[,,,P,,]
   V <- V[,,,P,,]

   ws.o  <- sqrt(U[,,,,obindex]^2  + V[,,,,obindex]^2)
   ws.sd <- sqrt(U[,,,,oeindex]^2  + V[,,,,oeindex]^2)
   u.f   <- U[,,,,obindex] - U[,,,,fcindex]
   v.f   <- V[,,,,obindex] - V[,,,,fcindex]
   ws.f  <- sqrt(u.f^2 + v.f^2)

   W <- V ## use to create array
   W[,,,,obindex] <- ws.o 
   W[,,,,fcindex] <- ws.o - ws.f 
   W[,,,,oeindex] <- ws.sd 

   # wind speed censor. Anything over 40 (m/s?) is unrealistic
   W[ W > 40] <- NA

} # close ws loop


#----------------------------------------------------------------------
# Determine the "id"s or R MARGINS for different tasks

## TJH FIXME - what is a reasonable default for group.dims
## Match the strings to the dimension indices - needed by 'apply'
## if group.dims == ("level","station") then group.id == 4,5

print('determining the group ids')

group.dims


group.id <- array(0,dim=length(group.dims))

for (i in 1:length(group.dims)){
   group.id[i] <- match(group.dims[i], org$dimstrings)
   if (is.na(group.id[i])){
      stop(sprintf("ERROR: cannot match group string %s",group.dims[i]))
   }
}

## seems more natural to have the ids in an ascending order, must
## permute the strings to match the ascending order.

bob        <- sort(group.id, index.return=T)
group.dims <- group.dims[bob$ix]
group.id   <- group.id[  bob$ix]

rm(bob)

print('Done determining the group ids.')

#----------------------------------------------------------------------
# bootstrap - from the set of analysis cycles (# of forecasts available)
# No need to bootstrap if there is only one analysis cycle ... true?
#----------------------------------------------------------------------
#### These need to be hard coded to reflect dimensions of ensemble and bias

OUT.bias <- OUT.thres <- OUT.rank <- OUT.bs <- OUT.attr <- OUT.crps <- list()

set.seed(12)

if (n.analysis < 10){
   warning(sprintf('Not enough analyses (%d) to be able to bootstrap',n.analysis))
   n.boot <- 1
}

for( g in 1:n.boot){

   print(sprintf('Trip %d of %d down bootstrap lane for %s',g,n.boot,filename))

   if(g == 1){ 
      analysis.id <- 1:n.analysis
   } else {
      analysis.id <- sample(1:n.analysis, n.analysis, replace = TRUE)
   }

   boot <- list(forecast   = org$forecast[,,,,analysis.id,drop=F],
                observed   = org$observed[,,,,analysis.id,drop=F],
                oberrvar   = org$oberrvar[,,,,analysis.id,drop=F],
                bias.dims  = bias.dims,
                dimstrings = org$dimstrings,
                dimlengths = org$dimlengths)

   #----------------------------------------------------------------------
   # The set of analysis.id's chosen for each bootstrap sample means 
   # a distinct bias should be removed.
   #----------------------------------------------------------------------

   if(length(bias.dims) > 0){

      print(sprintf('Trip %d of %d : removing bias.',g,n.boot))
      Wb <- bias.rm( boot )

      boot$forecast <- Wb$data  # bias free forecast
      OUT.bias[[g]] <- Wb$bias
      rm(Wb)

   } else {

      OUT.bias[[g]] <- NA

   }

   #----------------------------------------------------------------------
   # add observation error simulated with a normal distribution to forecasts
   #----------------------------------------------------------------------

   if(add.obs.error){

      # set.seed(11) ## useful if you want to replicate results.
      print(sprintf('Trip %d of %d : adding observation error',g,n.boot))

      oberrorSD <- sqrt(boot$oberrvar)
      xx        <- as.numeric(oberrorSD)
      yy        <- rnorm(n = length(xx), mean = 0, sd = xx)
      ZZ        <- array(yy, dim = dim(oberrorSD))

      boot$forecast <- boot$forecast + ZZ 	

      rm(ZZ,yy,xx,oberrorSD)
   }

   #----------------------------------------------------------------------
   # Calculate thresholds either [by level] or [by level, by site], or 
   # a [chosen threshold]. 
   # There are several ways to convert continous forecasts into binary events.
   # A threshold needs to be chosen.
   # A single threshold can be chosen for the entire domain.  
   # A quantile can be selected.
   # The quantile can be used to calculate different thresholds for level, stations, etc.

   # under the old storage paradigm
   #group.id <- c(4,5) ## when there is  a level coord : 5 corresponds to site index, 4 = level
   #group.id <- c(4)   ## when there is NO level coord : 4 corresponds to site index

   # TJH - fixme - from here on, it makes sense to just use 1 copy of the obs & obs error

   print(sprintf('Trip %d of %d : working on quantiles',g,n.boot))

   ## obs.binary and exc.o.thresh are logical arrays
   if(prob.method == "fixed"){

      obs.binary   <- boot$observed > fixed.thres
      exc.o.thres  <- boot$forecast > fixed.thres  

   } else {

      ## threshold by quantiles
      ## thresh   is the matrix of quantiles  (nlevel -x- nstation)

      if (length(quantiles) == 1) {

         thres        <- apply(boot$observed,group.id,quantile,quantiles[1],na.rm=TRUE)
         obs.binary   <- sweep(boot$observed,group.id,thres, FUN = ">")
         exc.o.thres  <- sweep(boot$forecast,group.id,thres, FUN = ">")

      } else {

         thres        <- apply(boot$observed,group.id,quantile,quantiles[1],na.rm=TRUE)
         obs.binary1  <- sweep(boot$observed,group.id,thres, FUN = ">")
         exc.o.thres1 <- sweep(boot$forecast,group.id,thres, FUN = ">")

         thres        <- apply(boot$observed,group.id,quantile,quantiles[2],na.rm=TRUE)
         obs.binary2  <- sweep(boot$observed,group.id,thres, FUN = "<")
         exc.o.thres2 <- sweep(boot$forecast,group.id,thres, FUN = "<")

         obs.binary   <- obs.binary1 & obs.binary2
         exc.o.thres  <- exc.o.thres1 & exc.o.thres2

      }  
   }

   # TJH DEBUG ... thres is not defined for a fixed threshold
   # storing threshold
   OUT.thres[[g]] <- thres	

   # convert to probabilisitic forecast and binary observation
   # sum(logical vector) is counting ... averaging over ensemble dimension

   # prob.frcst <- apply(exc.o.thres, c(1,2,4,5), sum) / length(ensemble)
   # prob.frcst <- apply(exc.o.thres, c(1,2,4),   sum) / length(ensemble)
   # 1,2,4,5 is 'time','date','level','station'
   # 1,2,4   is 'time','date','station'   (when there is no 'level')

   # This squeezes out the ensemble dimension
   prob.frcst <- apply(exc.o.thres, neworder[-id$ens], sum) / length(ensemble)

   # All the observation copies are the same, pick 1, match shape of prob.frcst
   # TJH - so why are we not subsetting 'observed' right off the bat?
   obs.binary <- array(obs.binary[,1,,,], dim=dim(prob.frcst))

   ### So by now, we have lost the ensemble dimension
   ### obs.binary ~ [n.lead  n.level  n.station  n.analysis]
   ### prob.frcst ~ [n.lead  n.level  n.station  n.analysis]

   dimstrings = org$dimstrings[-id$ens]

   #----------------------------------------------------------------------
   # TJH - dunno what this is
   #----------------------------------------------------------------------

   print(sprintf('Trip %d of %d : working on attr and BS',g,n.boot))

   if( run.attr ){

      nx = n.level * n.station * n.analysis

      ### temp1 has three 1D arrays for each forecast lead time 
      temp <- array(NA,dim=c(nx, 3*n.lead))

      for(z in 1:n.lead){ ## loop over forecast lead times
         temp[,3*(z-1)+1] = as.numeric(obs.binary[z,,,])
         temp[,3*(z-1)+2] = as.numeric(prob.frcst[z,,,])
         temp[,3*(z-1)+3] = rep(leads[z], nx)
      }

      OUT.attr[[g]] <- temp
   }	

   #----------------------------------------------------------------------k
   # BS
   #----------------------------------------------------------------------

   if(run.bs){

      temp        <- prob.frcst - obs.binary
      ff          <- function(x){mean(x^2, na.rm = TRUE)}
      OUT.bs[[g]] <- apply(temp, c(1), ff )

      rm(temp, nx)
   }

   #----------------------------------------------------------------------
   # crps
   #----------------------------------------------------------------------
  
   if(run.crps){

      print(sprintf('Trip %d of %d : working on CRPS',g,n.boot))

      mu    <- apply(boot$forecast, neworder[-id$ens], mean )
      sddat <- apply(boot$forecast, neworder[-id$ens], sd   )
      CRPS  <- array(NA,dim=n.lead)

      for(k in 1:n.lead){  ## loop in lead times
         ## all ensemble copies same, only need 1 for obs.mat
         mu.mat  <-       mu[k,,,]
         sd.mat  <-    sddat[k,,,]
         obs.mat <- array(boot$observed[k,1,,,], dim=dim(mu.mat))
         xx      <- crps(obs = as.numeric(obs.mat), 
                                  cbind(as.numeric(mu.mat), 
                                        as.numeric(sd.mat)) )
         CRPS[k] <- mean(xx$crps, na.rm = TRUE)
      }
      OUT.crps[[g]] <- CRPS

      rm(xx, obs.mat, sd.mat, mu.mat, CRPS, sddat, mu)

   }

   #----------------------------------------------------------------------
   # Calculate rank histogram.
   # TJH - I Think - I have no idea how this works ... or IF it works.
   # My test on the original data did not provide the result I was expecting.
   #----------------------------------------------------------------------
   ##                               n.lead n.anal n.ens n.site
   ## original dim(forecast) == [1]    6     32     10   3003
   ## original dim(temp.ff)  == [1]    6     32    3003   10
   ##                               n.lead n.anal n.site n.ens

   if(run.rank){

      print(sprintf('Trip %d of %d : rank histogram',g,n.boot))

      # permute so the last dimension is the ensemble.
      # TJH does it really matter if its the last (or the first?)

      myorder <- c(neworder[-id$ens], id$ens)
      temp.ff <- aperm(boot$forecast, perm = myorder) 
      temp.oo <- aperm(boot$observed, perm = myorder) 

      ## added dimension will wrap data (TJH - Explain to me).
      ## make a new array such that we can paste the observation onto
      ## the forecasts.  The rank function takes care of the rest
      nrhbins                <- n.ens + 1
      mydims                 <- dim(temp.ff)
      mydims[length(mydims)] <- nrhbins
      temp.ff2               <- array(temp.ff, dim = mydims)
      temp.ff2[,,,,nrhbins]  <- temp.oo[,,,,1]  ## add obs to last dimension

      f <- function(x){rank(x, ties.method = "random", na.last = "keep")[nrhbins] }

      OUT.rank[[g]] <- apply(temp.ff2, c(1,2,3), f)
   }
}  ## close boot loop

print('Done with bootstrap, saving output ...')

save.image(out.file)
}
