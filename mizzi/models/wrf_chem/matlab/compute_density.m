function rho = compute_density( mu, dnw, phi )
%% FUNCTION rho : Computes density of dry air, computed from the
%    hydrostatic relation for dry air: - dphi/deta = -mu alpha.
%
% Inputs:  
%	mu  = (full) dry hydrostatic surf. press (2d)
%	dnw = intervals between w levels, at mass pts
%	phi = (full) geopotential, at w pts
% Output:
%	rho  =  density, at mass pts
%
% See wrf subroutine calc_p_rho_phi.

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL: https://subversion.ucar.edu/DAReS/DART/branches/development/models/wrf/matlab/compute_density.m $
% $Id: compute_density.m 5074 2011-07-15 17:06:58Z thoar $
% $Revision: 5074 $
% $Date: 2011-07-15 11:06:58 -0600 (Fri, 15 Jul 2011) $

phi_eta = ( phi(2:end,:,:) - phi(1:end-1,:,:) ) ./ repmat( dnw(:), [1 size(mu)] );

rho = -repmat( reshape( mu, [1 size(mu)] ), [ length(dnw) 1 1 ] ) ./ phi_eta ;