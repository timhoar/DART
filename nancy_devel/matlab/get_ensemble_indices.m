function [ens_size, ens_indices] = get_ensemble_indices(fname)
%% DART:GET_ENSEMBLE_INDICES  returns the number of ensemble members in the file and their 'copy' indices. 
%
% Example:
% fname = 'Prior_Diag.nc';
% [ens_size, ens_indices] = get_ensemble_indices(fname);
%
% Example to return just the size ...
% [ens_size, ~] = get_ensemble_indices(fname);

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/mpas/matlab/get_ensemble_indices.m $
% $Id: get_ensemble_indices.m 5616 2012-03-22 22:42:39Z thoar $
% $Revision: 5616 $
% $Date: 2012-03-22 16:42:39 -0600 (Thu, 22 Mar 2012) $

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

ens_size = [];
ens_indices = [];

metastrings = nc_varget(fname,'CopyMetaData');
if(size(metastrings,2) == 1), metastrings = metastrings'; end
metadata = cellstr(metastrings);

% If the only copy is the true state, return without issuing the warning.

if strncmpi('true state',metadata,length('true state'))
   return
end

% see what we have ...

ens_indices = find(strncmpi('ensemble member',metadata,length('ensemble member')));

if (isempty(ens_indices))
   fprintf('WARNING: unable to find any valid ensemble members in %s\n', fname)
   disp('valid metadata strings are: ')
   for i = 1:length(metadata),
      fprintf('%s\n',metadata{i})
   end
   ens_size = [];
else
   ens_size = length(ens_indices);
end

