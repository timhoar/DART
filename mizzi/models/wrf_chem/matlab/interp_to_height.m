function var_interp = interp_to_height( var_in, heights, level)
%% interp_to_height Interpolates to a height level given heights
%
% var_interp = interp_to_height( var_in, heights, level)
%
% Interpolates var_in(:,:,:) to a height level given heights
% in heights(:,:,:).  Interpolation is linear in height. 
% Set var_interp to NaN where level is beneath (<) heights(1,:,:).

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL: https://subversion.ucar.edu/DAReS/DART/branches/development/models/wrf/matlab/interp_to_height.m $
% $Id: interp_to_height.m 5074 2011-07-15 17:06:58Z thoar $
% $Revision: 5074 $
% $Date: 2011-07-15 11:06:58 -0600 (Fri, 15 Jul 2011) $

[Nk Nj Ni] = size(heights); 
below      = zeros(Nj,Ni);
var_below  = zeros(Nj,Ni);
dz         = zeros(Nj,Ni);
dvar       = zeros(Nj,Ni);

% % at each horiz. location, find highest height beneath level
for ii = 1:Ni
for jj = 1:Nj

    % kk(jj,ii) = max( find( heights(:,jj,ii) - level < 0 ) ) ;
    kk = max( find( heights(:,jj,ii) - level < 0 ) ) ;
    % if isempty(kk); disp([ kk , jj, ii ] ); end
    if isempty(kk); kk = 1; end  % level is below surface here

       below(jj,ii) = heights(max(1 ,kk  ),jj,ii) ;  % height below level
          dz(jj,ii) = heights(min(Nk,kk+1),jj,ii) - heights(max(1 ,kk  ),jj,ii) ;
   var_below(jj,ii) = var_in( max(1 ,kk  ),jj,ii) ;
       dvar (jj,ii) = var_in( min(Nk,kk+1),jj,ii) -  var_in(max(1 ,kk  ),jj,ii) ;
end
end

var_interp =  ( dvar ./ dz ) .* ( level - below )  + var_below ;

var_interp( level < heights(1,:,:) ) = NaN ;
% level is beneath lowest mass level