function height = compute_height( phi, g )
%% compute_height
%
% Inputs:
%	phi = (full) geopotential, at w pts
%	g   = gravitational acceleration
% Output:
%	height = height, at mass pts
%

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL: https://subversion.ucar.edu/DAReS/DART/branches/development/models/wrf/matlab/compute_height.m $
% $Id: compute_height.m 5074 2011-07-15 17:06:58Z thoar $
% $Revision: 5074 $
% $Date: 2011-07-15 11:06:58 -0600 (Fri, 15 Jul 2011) $

height = ( phi(2:end,:,:) + phi(1:end-1,:,:) ) ./ (2*g) ;
