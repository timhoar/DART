function temp = compute_temperature( pres, theta, Cp, Rd, p0 )
%% FUNCTION compute_temperature - Computes temperature from potential temperature.
%
% Inputs:  
%	pres     = pressure, at mass pts
%       theta    = potential temperature, at mass pts 
%	Cp,Rd,p0 = c_p, dry gas constants, reference pressure
% Output:
%	temp     = temperature, at mass pts

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL: https://subversion.ucar.edu/DAReS/DART/branches/development/models/wrf/matlab/compute_temperature.m $
% $Id: compute_temperature.m 5074 2011-07-15 17:06:58Z thoar $
% $Revision: 5074 $
% $Date: 2011-07-15 11:06:58 -0600 (Fri, 15 Jul 2011) $

kappa = Rd / Cp ;

temp = theta .* (pres ./ p0).^kappa ;
