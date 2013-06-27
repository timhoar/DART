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
# Creates new ncdf files using attributes in the globabl environment.
# This is a little sloppy.
#
#======================================================================

f.out <- function() {
	
a <- files[j]
	
out.name <- paste(out.dir, var.char, substr(a,1,nchar(a)-3), ".out.nc", sep = "")
	
date.dim             <- ex.nc$dim$date
date.dim$longname    <- "date"

station.dim          <- ex.nc$dim$station
station.dim$longname <- "station" 
station.dim$value    <- station  ### correct for length of surface stations
station.dim$len      <- length(station)  ### correct for length of surface stations

level.dim            <- ex.nc$dim$level
level.dim$longname   <- "Level" 

ensemble.dim          <- ex.nc$dim$ensemble
ensemble.dim$longname <- "ensemble id"

time.dim              <- ex.nc$dim$time
time.dim$longname     <- "Lead Time hours"
time.dim$vals         <- seq(0,60, 12) ### does not need to be hard coded if files are correct.

TMP.list <- list(time.dim, date.dim, ensemble.dim, level.dim, station.dim)

mv     <- -888888

bias   <- var.def.ncdf("bias", units = "double", 
                       dim = list(time.dim, ensemble.dim, level.dim, station.dim),  
                       missval = mv, 
                       longname = "Bias calcuated at each point")

mse    <- var.def.ncdf("mse", units = "double", 
                       dim = list(time.dim, ensemble.dim, level.dim, station.dim),
                       missval = mv,
                       longname = "mean squared error")

n      <- var.def.ncdf("n", units = "counts", 
                       dim = list(time.dim, ensemble.dim, level.dim, station.dim),  
                       missval = mv, longname = "Number of data points at each location")

if(prob.method == "quantile"){

   att.put.ncdf( new, 0, "prob.method", "quantile")
   p.thres <- var.def.ncdf("p.thres", units = v.units, 
                       dim = TMP.list[group.id],  
                       missval = mv, 
                       longname = "Thresold used to create probabilistic forecast")

   var.add.ncdf(new, v = p.thres)
   close.ncdf(new)

   new    <- open.ncdf(out.name, write = TRUE)

   put.var.ncdf(new, varid = p.thres, vals = thres )

}
	
if(prob.method == "fixed"){
   att.put.ncdf( new, 0, "prob.method", "fixed")
   p.thres <- var.def.ncdf("p.thres", units = v.units,
                       dim = 1,
                       missval = mv,
                       longname = "Thresold used to create probabilistic forecast")
}	

prob.frcst.ncdf <- var.def.ncdf("prob.forecast", units = "Probability", 
                       dim = list( time.dim, ensemble.dim, level.dim, station.dim),  
                       missval = mv, 
                       longname = "Probability of exceeding threshold")

prob.obs.ncdf   <- var.def.ncdf("prob.obs", units = "Binary",
                       dim = list( time.dim, ensemble.dim, level.dim, station.dim),
                       missval = mv,
                       longname = "Outcome of event")

#### fill variable

new <- create.ncdf( out.name, list(bias, mse, n, p.thres) )

put.var.ncdf(new, varid = prob.frcst.ncdf, vals = prob.frcst)
put.var.ncdf(new, varid = prob.frcst.obs,  vals = obs.binary)
put.var.ncdf(new, varid = p.thres,         vals = thres)
put.var.ncdf(new, varid = bias,            vals = Wb$bias)
put.var.ncdf(new, varid = mse,             vals = Wb$mse)
put.var.ncdf(new, varid = n,               vals = Wb$n)

att.put.ncdf( new, 0, "source", NMS)
att.put.ncdf(new, p.thres, "prob.method", "quantile")
close.ncdf(new)

}

