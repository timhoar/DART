function var_interp = interp_to_pressure( var_in, pressure, p_level)
% Interpolates var_in(:,:,:) to a pressure level p_level given pressures
% in pressures(:,:,:).  Interpolation is linear in log pressure. 
% Set var_interp to NaN where p_level is beneath (>) pressure(1,:,:).

% Data Assimilation Research Testbed -- DART
% Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

% <next three lines automatically updated by CVS, do not edit>
% $Id$
% $Source$
% $Name$

 [Nk Nj Ni] = size(pressure) ; 

 below = zeros(Nj,Ni) ;  var_below = zeros(Nj,Ni) ;
 dlogp = zeros(Nj,Ni) ; dvar = zeros(Nj,Ni) ;

 % at each horiz. location, find highest pressure level beneath p_level
 for ii = 1:Ni; for jj = 1:Nj; 
    % kk(jj,ii) = max( find( pressure(:,jj,ii) - p_level > 0 ) ) ;
    kk = max( find( pressure(:,jj,ii) - p_level > 0 ) ) ;
    % if isempty(kk); disp([ kk , jj, ii ] ); end
    if isempty(kk); kk = 1; end  % p_level is below surface here

    below(jj,ii) = log( pressure(max(1 ,kk  ),jj,ii) ) ;  % pressure level below p_level
    dlogp(jj,ii) = log( pressure(min(Nk,kk+1),jj,ii) )    ...
                    - log( pressure(max(1 ,kk  ),jj,ii) ) ;
    var_below(jj,ii) = var_in(max(1 ,kk  ),jj,ii) ;
    dvar (jj,ii) = var_in(min(Nk,kk+1),jj,ii) - var_in(max(1 ,kk  ),jj,ii) ;
 end; end

 var_interp =  ( dvar ./ dlogp ) .* ( log(p_level) - below )  + var_below ;

 var_interp( p_level > pressure(1,:,:) ) = NaN ;
   % p_level is beneath surface