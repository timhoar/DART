function slab = get_hyperslab(varargin)
%% DART:get_hyperslab gets an arbitrarily-shaped variable from a netCDF diagnostic file.
% Retrieves a particular copy of a state vector from a file whose
% full or relative path is specified in the file argument.

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/mpas/matlab/get_hyperslab.m $
% $Id: get_hyperslab.m 5616 2012-03-22 22:42:39Z thoar $
% $Revision: 5616 $
% $Date: 2012-03-22 16:42:39 -0600 (Thu, 22 Mar 2012) $

for i = 1:2:nargin,
   eval(sprintf('pinfo.%s = varargin{i+1};',varargin{i}))
end

if ( exist(pinfo.fname,'file') ~= 2 ), error('%s does not exist.',pinfo.fname); end

[start, count] = GetNCindices(pinfo,'fname',pinfo.varname);
slab           = nc_varget(pinfo.fname, pinfo.varname, start, count);

if (sum(isfinite(slab(:))) == 0)
   pinfo
   error('%s %s has all missing values ... exiting.', pinfo.fname, pinfo.varname)
end

