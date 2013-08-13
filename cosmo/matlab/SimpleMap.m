function h = SimpleMap(fname, varname, time, copystring, level)
%% SimpleMap plots a horizontal slice from a DART netcdf file.
%
% Since WRF does not have tidy coordinate variables for lat/lon,
% this script will not directly support WRF/DART netCDF files.
%
% EXAMPLE:
%
% fname   = 'PriorDiag.nc';
% varname = 'T';      % 'PS','T','US','VS','Q','CLDLIQ','CLDICE' ...
% time    = 1;        % index into time array
% copy    = 'ensemble mean';  % anything from ncdump -v CopyMetaData
% level   = 1;        % index into level array
% h = SimpleMap(fname,varname,time,copy,level);
%
% EXAMPLE:
%
% fname   = '/project/dart/raeder/CSL/Groups2/PriorDiag.nc';
% varname = 'T';      % 'PS','T','US','VS','Q','CLDLIQ','CLDICE' ...
% time    = 1;        % index into time array
% copy    = 'ensemble mean';  % anything from ncdump -v CopyMetaData
% level   = 1;        % index into level array
% h = SimpleMap(fname,varname,time,copy,level);

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist(fname,'file') ~= 2), error('%s does not exist.',fname); end

% Get some information from the netCDF file.

plotdat = netCDFinfo(fname, varname);

% Figure out the starting date

timebase   = sscanf(plotdat.timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin = datenum(timebase(1),timebase(2),timebase(3));
timearr    = nc_varget(fname,'time');
timestring = datestr(timearr(1) + timeorigin);

% create the hyperslab indices ... 

copyindex = get_copy_index(fname,copystring);

myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
myinfo.levelindex = level;
myinfo.timeindex  = time;
[start, count]    = GetNCindices(myinfo,'diagn',varname);

datmat = nc_varget(fname, varname, start, count);

% Create the plot.
% The figure window is partitioned into two parts ... one has 
% the actual graphic ... the other (the bottom part) is a text
% table of all the attributes of interest for the graphic.

h0 = subplot('position',[0.1 0.3 0.8 0.6]); set(h0,'FontSize',14)
h = imagesc(plotdat.lons,plotdat.lats,datmat);
set(gca,'YDir','normal');
axis image;
worldmap;
h1 = colorbar;
set(h1,'FontSize',14)
set(get(h1,'Ylabel'),'String',plotdat.varunits,'FontSize',14)
title(sprintf('%s ''%s''',plotdat.model,varname))

t0 = subplot('position',[0.1 0.05 0.8 0.25]);
t1 = plot([0 1],[0 1],'.','Visible','off');
axis off

dimstring = [];
for i = 1:length(plotdat.dimnames)
   dimstring = sprintf('%s %s',dimstring,plotdat.dimnames{i});
end

text(0.0, 0.65, sprintf('%s %s',plotdat.model,fname),'Interpreter','none');
text(0.0, 0.50, timestring);
text(0.0, 0.35, sprintf('%s ''%s'' %s [%s]', ...
               copystring, varname, plotdat.varlname, dimstring));
text(0.0, 0.20, sprintf('level %.2f (index %d) at timestep %d', ...
                         plotdat.levels(level), level, time));

annotate(plotdat)



% helper function(s) follow ...

function plotdat = netCDFinfo(fname,varname)

plotdat.model      = nc_attget(fname, nc_global, 'model');
plotdat.varlname   = nc_attget(fname, varname,   'long_name');
plotdat.varunits   = nc_attget(fname, varname,   'units');
plotdat.timeunits  = nc_attget(fname,'time',     'units');
plotdat.calendar   = nc_attget(fname,'time',     'calendar');

varinfo            = nc_getvarinfo(fname,varname);
plotdat.dimnames   = varinfo.Dimension;

for i = 1:length(varinfo.Dimension)

   diminfo = nc_getdiminfo(fname,varinfo.Dimension{i});
   dimname = diminfo.Name;

   % the coordinate variable is the same name as the dimension
   % some grids have multiple vertical levels so there is no one
   % predictable coordinate variable name. By convention, these
   % variables should be tagged with a 'cartesian_axis' attribute.
   % ditto for lat, lon ... (on staggered grids, for example)
   % So the XG coordinate dimension has 'cartesian_axis = X',
   % for example.

   [len, status, value] = is_dimension_cartesian(fname, diminfo.Name);

   if (status > 0)
      dimname = value;
   else
      % Then there is no 'cartesian_axis' attribute and the best we can
      % hope for is a standard dimension name [time,copy,lat,lon,lev]
   end

   switch lower(dimname) % loop over all likely coordinate variables
      case {'lev','z'}
         plotdat.levdimname = diminfo.Name;
         plotdat.levels     = nc_varget(fname,diminfo.Name);
         plotdat.levunits   = nc_attget(fname,diminfo.Name,'units');
         plotdat.levlname   = nc_attget(fname,diminfo.Name,'long_name');
      case {'lat','y','tmpj'}
         plotdat.latdimname = diminfo.Name;
         plotdat.lats       = nc_varget(fname,diminfo.Name);
         plotdat.latunits   = nc_attget(fname,diminfo.Name,'units');
      case {'lon','x','tmpi'}
         plotdat.londimname = diminfo.Name;
         plotdat.lons       = nc_varget(fname,diminfo.Name);
         plotdat.lonunits   = nc_attget(fname,diminfo.Name,'units');
   end
end



function [len, status, value] = is_dimension_cartesian(fname,dimname)

status   = 0;
len      = 0;
value    = [];
Cvarinfo = nc_getvarinfo(fname, dimname);

for j = 1:length(Cvarinfo.Attribute);
   attribute = Cvarinfo.Attribute(j);
   switch lower(attribute.Name)
      case{'cartesian_axis'}
         status = 1;
         len    = Cvarinfo.Size;
         value  = attribute.Value;
         break
      otherwise
   end
end



function annotate(plotdat)

nlevels = length(plotdat.levels);

if (nlevels == 1) 
   return
end

stride = 4;

dx = 1.0/ceil(nlevels/stride);

text(0, 0.05,'level index','FontSize',8)
text(0,-0.05,'level value','FontSize',8)
text(0,-0.15,plotdat.levlname,'FontSize',8)

itrip = 1;
for i = 1:stride:nlevels
   h1 = text(dx*itrip, 0.05,sprintf('%d',                  i ));
   h2 = text(dx*itrip,-0.05,sprintf('%9.2f',plotdat.levels(i)));
   itrip = itrip + 1;
   set(h1,'FontSize',8,'HorizontalAlignment','right')
   set(h2,'FontSize',8,'HorizontalAlignment','right')
end