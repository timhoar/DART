function h = map_wrf(fname, varname, levelindx, timeindx, copystring )
%% map_wrf creates an image of a WRF field without using the mapping toolbox.
%
%% Example using a DART-style diagnostic file, i.e.:
%
% fname      = '/glade/proj3/image/hliu/200812new/cwb_icbc/12_01/Prior_Diag.nc';
% varname    = 'U_d01';
% copystring = 'ensemble mean';
% map_wrf(fname, varname, levelindx, timeindx, copystring);
% worldmap;
% axis off;
%
%
%% Example using a WRF netCDF file:
%
% fname      = '/glade/proj3/image/hliu/ICBC_from_cwb/wrfinput_d01';
% levelindx  = 10;
% timeindx   = 1;
% map_wrf(fname, 'U', levelindx, timeindx);
% worldmap;   % superimpose some low-res coastlines 
% axis off;
%
%% layer on the locations of some observations:
%
% obsfile = 'obs_epoch_001.nc';
% ObsTypeString = 'RADIOSONDE_U_WIND_COMPONENT';
% region        = [0 360 -90 90 -Inf Inf];
% CopyString    = 'observation';
% QCString      = 'DART quality control';
% verbose       = 1;   % anything > 0 == 'true'
% obs = read_obs_netcdf(obsfile, ObsTypeString, region, CopyString, QCString, verbose);
% hold on;
% plot(obs.lons, obs.lats, 'kd', 'MarkerFaceColor','k')

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% $Id$

copystring = 'no copy';
levelindex = 1;
timeindex  = 1;

for iarg = 1:2:size(varargin,2)
   mystring = char(varargin{iarg});
% debug fprintf('argument %d is %s\n',iarg,mystring)
   switch lower(mystring(1:3))
   case 'lev'
      levelindex = round(varargin{iarg+1});
   case 'cop'
      copystring =  char(varargin{iarg+1});
   case 'tim'
      timeindex  = round(varargin{iarg+1});
   otherwise
      varargin{iarg}
      error('unknown argument')
   end
end

% Check to ensure the file and the desired variable exist

if (exist(fname,'file') ~= 2) 
   error('%s does not exist',fname)
end

varexist(fname,{varname})

%% Get the dimension information.

xy          = GetWRFlatlon( fname, varname);
levels      = GetWRFlevels( fname, varname);
times       = GetWRFtimes(  fname, varname);
timestrings = datestr(times,31);

%% Set the projection parameters
%  float MAP_PROJ(domain) ;
%        MAP_PROJ:long_name = "domain map projection" ;
%        MAP_PROJ:units = "0=none, 1=Lambert, 2=polar, 3=Mercator, 5=Cylindrical, 6=Cassini" ;

if     ( mystruct.mp(domain) == 0 ) 

   error('unsupported projection code %d',mystruct.mp(domain))

elseif ( mystruct.mp(domain) == 1 ) % a conic projection
    % projection attributes
    mystruct.map_proj    = 'lambert';
    mystruct.maplatlimit = [min(mystruct.yvar(:)) max(mystruct.yvar(:))];
    mystruct.maplonlimit = [min(mystruct.xvar(:)) max(mystruct.xvar(:))];
    mystruct.clongitude  = mystruct.stand_lon(domain); % some central longitude
    mystruct.parallels   = [mystruct.stdlat1(domain), mystruct.stdlat2(domain)];
    mystruct.gridargs    = {'box','off','xaxis','top','xtick',12};

elseif ( mystruct.mp(1) == 2 ) 
    % projection attributes
    mystruct.map_proj     = 'ups';
    mystruct.maplatlimit  = [minlat maxlat];
    mystruct.maplonlimit  = [minlon maxlon];
    mystruct.radius       = double(ceil(abs(diff(mystruct.maplatlimit))));
    if (absmin == minlat)
        mystruct.zone   = 'north';
        mystruct.origin = 90;
    else
        mystruct.zone   = 'south';
        mystruct.origin = -90;
    end
    mystruct.gridargs    = {'box','off','xaxis','top','xtick',12};

else

   error('unknown projection code %d',mystruct.mp(domain))

end



mystruct.timevarname  = 'time';
mystruct.time         = nc_varget(fname, 'time');
mystruct.timeunits    = nc_attget(fname, 'time', 'units');
mystruct.calendar     = nc_attget(fname, 'time', 'calendar');
mystruct.ntimes       = length(mystruct.time);     %   YYYY   MM   DD    HH    MI     S
timebase              = sscanf(mystruct.timeunits,'%*s%*s%d%*c%d%*c%d%*c%2d%*c%2d%*c%2d');
mystruct.timeorigin   = datenum(timebase(1),timebase(2),timebase(3), ...   % YYYY MM DD
                                timebase(4),timebase(5),timebase(6));      % H, MI, S
mystruct.timestrings  = datestr(mystruct.time + mystruct.timeorigin, 0);

%% Get the metadata for the variable 

varinfo = nc_getvarinfo(fname, varname);

for i = 1:length(varinfo.Attribute)

   attname = varinfo.Attribute(i).Name;
   switch lower(attname)
      case 'long_name'
         mystruct.varlongname = varinfo.Attribute(i).Value;
      case 'units'
         mystruct.varunits = varinfo.Attribute(i).Value;
   end
end

copydim = find(strncmp('copy',varinfo.Dimension,length('copy')));

if (isempty(copydim))
   copyindex = NaN;
else
   copyindex = get_copy_index(fname,copystring);
end

for itime = timeindx

   plot_title = {fname, ...
         sprintf('levelindex %d %s %s',levelindx,varname,timestrings(itime,:))};

   %% Determine the hyperslab indexing
   
   myinfo.diagn_file = fname;
   if (isfinite(copyindex)), myinfo.copyindex  = copyindex; end
   myinfo.levelindex = levelindx;
   myinfo.timeindex  = timeindx;
   [start, count]    = GetNCindices(myinfo, 'diagn', varname);
   
% Extract field
   
   datmat = double(nc_varget(fname, varname, start, count));
   
   clf; 
   h = pcolor(xy.lonmat, xy.latmat, datmat);
   set(h,'LineStyle','none');
   
   title(plot_title,'Interpreter','none')
   h2 = colorbar('vert');
   set(get(h2,'YLabel'),'String',varinfo.units)

bob = get_var_grid(fname, varname);

if (isfield(bob,'xvarname')), mystruct.xvarname = bob.xvarname; end
if (isfield(bob,'xvar'    )), mystruct.xvar     = bob.xvar;     end
if (isfield(bob,'yvarname')), mystruct.yvarname = bob.yvarname; end
if (isfield(bob,'yvar'    )), mystruct.yvar     = bob.yvar;     end
if (isfield(bob,'zvar'    ))
   mystruct.zvarname  = bob.zvarname;
   mystruct.zvar      = bob.zvar;
   mystruct.level     = mystruct.zvar(mystruct.levelindex);
end
if ( isfield(bob,'zvarunits'   )), mystruct.zvarunits    = bob.zvarunits;    end
if ( isfield(bob,'zvarlongname')), mystruct.zvarlongname = bob.zvarlongname; end

%% Set the latitude/longitude limits

absmin = min(abs(mystruct.yvar(:)));
minlat = min(mystruct.yvar(:));
maxlat = max(mystruct.yvar(:));
minlon = min(mystruct.xvar(:));
maxlon = max(mystruct.xvar(:));

%% Read map metadata
% WRF name                                          Matlab name
% Cylindrical Equidistant Lat/Lon (code = 0)     eqdcylin (cylindrical
% Lambert Conformal               (code = 1)     lambert  (conic)
% Polar Stereographic             (code = 2)     ups      (azimuthal)
% Mercator                        (code = 3)     mercator (cylindrical)
% Gaussian                        (code = 4)
% Lat/Lon                         (code = 6)
% Rotated Lat/Lon                 (code = 203)

if (mystruct.mp == 0)
    % projection attributes
    mystruct.map_proj     = 'eqdcylin';
    mystruct.maplatlimit  = [-90 90];
    mystruct.maplonlimit  = [-180 180];
    mystruct.mapparallels = 30;
    mystruct.nparallels   = 1;
    mystruct.origin       = [0 0 0];
    % labelling attributes
    mystruct.mlinelocation  = 30;
    mystruct.mlinevisible   =  'on';
    mystruct.plinelocation  =  15;
    mystruct.plinevisible   = 'on';
    mystruct.labelformat    = 'compass';
    mystruct.labelrotation  = 'off';
    mystruct.labelunits     = 'degrees';
    mystruct.meridianlabel  = 'off';
    mystruct.mlabellocation = 30;
    mystruct.mlabelparallel = 90;
    mystruct.parallellabel  = 'off';
    mystruct.plabellocation = 15;
    mystruct.plabelmeridian = -180;

for i = 1:nvars
   gotone(i) = nc_isvar(filename,varnames{i});
   if ( ~ gotone(i) )
      warning('\n%s is not a variable in %s\n',varnames{i},filename)
   end
end

    mystruct.gridargs = {'box','off','xaxis','top','xtick',12};
%   m_grid('linewi',1,'tickdir','in','box','off','ytick',[-45 -60 -75], ...
%          'xaxis','top','xtick',12);
    
elseif (mystruct.mp == 3)
    mystruct.map_proj = 'mercator';
    mystruct.gridargs = {'box','off','xaxis','top','xtick',12};

else
    error('Unknown map projection ... %d',mystruct.mp)
end



function xy = GetWRFlatlon(fname, varname);
%% Each of the WRF variables has a 'coordinate' attribute signifying which
% of the 6 possible lat/lon variables is appropriate.

coordinates{1} = sscanf(nc_attget(fname,varname,'coordinates'),'%s %*s');
coordinates{2} = sscanf(nc_attget(fname,varname,'coordinates'),'%*s %s');

latcoord       = find(strncmp('XLAT', coordinates, length('XLAT')) > 0);
loncoord       = find(strncmp('XLON', coordinates, length('XLON')) > 0);
xy.latmat      = nc_varget(fname, coordinates{latcoord});
xy.lonmat      = nc_varget(fname, coordinates{loncoord});
xy.latunits    = nc_attget(fname, coordinates{latcoord},'units');
xy.lonunits    = nc_attget(fname, coordinates{latcoord},'units');

inds = (xy.lonmat < 0); % Convert to 0,360 to minimize dateline probs.
xy.lonmat(inds) = xy.lonmat(inds) + 360.0;



function levels = GetWRFlevels(fname, varname);
%% Return the appropriate vertical indices.

varinfo = nc_getvarinfo(fname,varname);
lvlind  = find(strncmp('bottom_top',varinfo.Dimension,length('bottom_top')));
nlevels = varinfo.Size(lvlind);
levels  = 1:nlevels;



function times = GetWRFtimes(fname, varname);
%% Return the appropriate time coordinates.
%  DART files use 'time' as a real array
%  WRF  files use 'Time' as a character string

varinfo = nc_getvarinfo(fname,varname);
timeind = find(strncmp('time',varinfo.Dimension,length('time')));

if (isempty(timeind)) % WRF flavor Time variable.
    timeind = find(strncmp('Time',varinfo.Dimension,length('Time')));
    if (isempty(timeind)), error('%s has no time information. Dying.'); end
    timestr  = nc_varget(fname, 'Times');  % Note the plural ... seeth.
    timebase = sscanf(timestr,'%d%*c%d%*c%d%*c%d%*c%d%*c%d'); % YYYY MM DD HH MM SS
    times    = datenum(timebase');

else
    times       = nc_varget(fname, 'time');
    calendar    = nc_attget(fname, 'time','calendar');
    timeunits   = nc_attget(fname, 'time','units');
    timebase    = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
    timeorigin  = datenum(timebase(1),timebase(2),timebase(3));
    times       = times + timeorigin;
 
end

% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$
