function pinfo = GetTIEGCMInfo(pstruct,fname,routine);
%% GetTIEGCMInfo   prepares a structure of information needed by the subsequent "routine"
%                The information is gathered via rudimentary "input" routines.
%
% pinfo = GetTIEGCMInfo(pstruct,routine);
%
% pstruct   structure containing the names of the truth_file and the diagn_file of the DART netcdf file
% routine   name of subsequent plot routine.

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

pinfo  = pstruct;
model  = nc_attget(fname,nc_global,'model');

if strcmpi(model,'TIEGCM') ~= 1
   error('Not so fast, this is not a TIEGCM model.')
end

%% Get the domain-independent information.

varexist(fname, {'copy','time'})

copy   = nc_varget(fname,'copy');
times  = nc_varget(fname,'time');

% Coordinate between time types and dates

timeunits  = nc_attget(fname,'time','units');
timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin = datenum(timebase(1),timebase(2),timebase(3));
dates      = times + timeorigin;

ilevel = nc_varget(fname,'ilev');    % interfaces
levels = nc_varget(fname, 'lev');    % midpoints
lon    = nc_varget(fname, 'lon');
lat    = nc_varget(fname, 'lat');


inds = find(lon >= 180);
lon(inds) = lon(inds) - 360.0;

prognostic_vars = get_DARTvars(fname);
num_vars = length(prognostic_vars);

switch lower(deblank(routine))

   case {'plotbins','plotenserrspread','plotensmeantimeseries','plotenstimeseries'}

      pgvar           = GetVar(prognostic_vars);  % Determine prognostic variable
      [level, lvlind] = GetLevel(pgvar,levels);   % Determine level and index
      [lat  , latind] = GetLatitude( pgvar,lat);
      [lon  , lonind] = GetLongitude(pgvar,lon);

      pinfo.model      = model;
      pinfo.times      = dates;
      pinfo.var        = pgvar;
      pinfo.level      = level;
      pinfo.levelindex = lvlind;
      pinfo.longitude  = lon;
      pinfo.lonindex   = lonind;
      pinfo.latitude   = lat;
      pinfo.latindex   = latind;


   case 'plotcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVar(prognostic_vars);
      [base_time, base_tmeind] = GetTime(dates);
      [base_lvl,  base_lvlind] = GetLevel(    base_var,levels);
      [base_lat,  base_latind] = GetLatitude( base_var,lat);
      [base_lon,  base_lonind] = GetLongitude(base_var,lon);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(prognostic_vars,          base_var);
      [comp_lvl, comp_lvlind] = GetLevel(    comp_var,levels,    base_lvlind);

      pinfo.model       = model;
      pinfo.times       = dates;
      pinfo.base_var    = base_var;
      pinfo.comp_var    = comp_var;
      pinfo.base_time   = base_time;
      pinfo.base_tmeind = base_tmeind;
      pinfo.base_lvl    = base_lvl;
      pinfo.base_lvlind = base_lvlind;
      pinfo.base_lat    = base_lat;
      pinfo.base_latind = base_latind;
      pinfo.base_lon    = base_lon;
      pinfo.base_lonind = base_lonind;
      pinfo.comp_lvl    = comp_lvl;
      pinfo.comp_lvlind = comp_lvlind;


   case 'plotvarvarcorrel'

      disp('Getting information for the ''base'' variable.')
       base_var                = GetVar(prognostic_vars);
      [base_time, base_tmeind] = GetTime(dates);
      [base_lvl , base_lvlind] = GetLevel(    base_var,levels);
      [base_lat , base_latind] = GetLatitude( base_var,lat);
      [base_lon , base_lonind] = GetLongitude(base_var,lon);

      disp('Getting information for the ''comparison'' variable.')
       comp_var               = GetVar(prognostic_vars,          base_var);
      [comp_lvl, comp_lvlind] = GetLevel(    comp_var,levels,    base_lvlind);
      [comp_lat, comp_latind] = GetLatitude( comp_var,lat,base_lat);
      [comp_lon, comp_lonind] = GetLongitude(comp_var,lon,base_lon);

      pinfo.model       = model;
      pinfo.times       = dates;
      pinfo.base_var    = base_var;
      pinfo.comp_var    = comp_var;
      pinfo.base_time   = base_time;
      pinfo.base_tmeind = base_tmeind;
      pinfo.base_lvl    = base_lvl;
      pinfo.base_lvlind = base_lvlind;
      pinfo.base_lat    = base_lat;
      pinfo.base_latind = base_latind;
      pinfo.base_lon    = base_lon;
      pinfo.base_lonind = base_lonind;
      pinfo.comp_lvl    = comp_lvl;
      pinfo.comp_lvlind = comp_lvlind;
      pinfo.comp_lat    = comp_lat;
      pinfo.comp_latind = comp_latind;
      pinfo.comp_lon    = comp_lon;
      pinfo.comp_lonind = comp_lonind;


   case 'plotsawtooth'

       pgvar          = GetVar(prognostic_vars);
      [level, lvlind] = GetLevel(    pgvar,levels);
      [  lat, latind] = GetLatitude( pgvar,lat);
      [  lon, lonind] = GetLongitude(pgvar,lon);
 %    [  lon, lonind] = GetCopies(pgvar,xxx);

      pinfo.model          = model;
      pinfo.times          = dates;
      pinfo.var_names      = pgvar;
      pinfo.truth_file     = [];
      pinfo.prior_file     = pstruct.prior_file;
      pinfo.posterior_file = pstruct.posterior_file;
      pinfo.level          = level;
      pinfo.levelindex     = lvlind;
      pinfo.latitude       = lat;
      pinfo.latindex       = latind;
      pinfo.longitude      = lon;
      pinfo.lonindex       = lonind;
      pinfo.copies         = 0;
      pinfo.copyindices    = [];

      if ( exist(pstruct.truth_file,'file') )
         pinfo.truth_file = pstruct.truth_file;
      end


   case 'plotphasespace'

      disp('Getting information for the ''X'' variable.')
       var1                   = GetVar(prognostic_vars);
      [var1_lvl, var1_lvlind] = GetLevel(    var1, levels);
      [var1_lat, var1_latind] = GetLatitude( var1, lat );
      [var1_lon, var1_lonind] = GetLongitude(var1, lon );

      disp('Getting information for the ''Y'' variable.')
       var2                   = GetVar(prognostic_vars,        var1    );
      [var2_lvl, var2_lvlind] = GetLevel(    var2, levels,     var1_lvlind);
      [var2_lat, var2_latind] = GetLatitude( var2, lat, var1_lat);
      [var2_lon, var2_lonind] = GetLongitude(var2, lon, var1_lon);

      disp('Getting information for the ''Z'' variable.')
       var3                   = GetVar(prognostic_vars,        var1    );
      [var3_lvl, var3_lvlind] = GetLevel(    var3, levels,     var1_lvlind);
      [var3_lat, var3_latind] = GetLatitude( var3, lat, var1_lat);
      [var3_lon, var3_lonind] = GetLongitude(var3, lon, var1_lon);

      % query for ensemble member
      s1 = input('Input ensemble member metadata STRING. <cr> for ''true state''  ','s');
      if isempty(s1), ens_mem = 'true state'; else ens_mem = s1; end

      % query for line type
      s1 = input('Input line type string. <cr> for ''k-''  ','s');
      if isempty(s1), ltype = 'k-'; else ltype = s1; end

      pinfo.model       = model;
      pinfo.times       = dates;
      pinfo.var1name    = var1;
      pinfo.var2name    = var2;
      pinfo.var3name    = var3;
      pinfo.var1_lvl    = var1_lvl;
      pinfo.var1_lvlind = var1_lvlind;
      pinfo.var1_lat    = var1_lat;
      pinfo.var1_latind = var1_latind;
      pinfo.var1_lon    = var1_lon;
      pinfo.var1_lonind = var1_lonind;
      pinfo.var2_lvl    = var2_lvl;
      pinfo.var2_lvlind = var2_lvlind;
      pinfo.var2_lat    = var2_lat;
      pinfo.var2_latind = var2_latind;
      pinfo.var2_lon    = var2_lon;
      pinfo.var2_lonind = var2_lonind;
      pinfo.var3_lvl    = var3_lvl;
      pinfo.var3_lvlind = var3_lvlind;
      pinfo.var3_lat    = var3_lat;
      pinfo.var3_latind = var3_latind;
      pinfo.var3_lon    = var3_lon;
      pinfo.var3_lonind = var3_lonind;
      pinfo.ens_mem     = ens_mem;
      pinfo.ltype       = ltype;

   otherwise

end


function pgvar = GetVar(prognostic_vars, defvar)
%----------------------------------------------------------------------
if (nargin == 2), pgvar = defvar; else pgvar = prognostic_vars{1}; end

str = sprintf(' %s ',prognostic_vars{1});
for i = 2:length(prognostic_vars),
   str = sprintf(' %s %s ',str,prognostic_vars{i});
end
fprintf('Default variable is ''%s'', if this is OK, <cr>;\n',pgvar)
fprintf('If not, please enter one of: %s\n',str)
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), pgvar = strtrim(varstring); end
inds        = strfind(pgvar,',');
pgvar(inds) = '';



function [time, timeind] = GetTime(times, deftime)
%----------------------------------------------------------------------
% Query for the time of interest.

ntimes = length(times);

% Determine a sensible default.
if (nargin == 2),
   time = deftime;
   tindex = find(times == deftime);
else
   if (ntimes < 2)
      tindex = round(ntimes/2);
   else
      tindex = 1;
   end
   time = times(tindex);
end

fprintf('Default time is %s (index %d), if this is OK, <cr>;\n',datestr(time),tindex)
fprintf('If not, enter an index between %d and %d \n',1,ntimes)
fprintf('Pertaining to %s and %s \n',datestr(times(1)),datestr(times(ntimes)))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), tindex = str2num(varstring); end

timeinds = 1:ntimes;
d        = abs(tindex - timeinds); % crude distance
ind      = find(min(d) == d);      % multiple minima possible
timeind  = ind(1);                 % use the first one
time     = times(timeind);



function [level, lvlind] = GetLevel(pgvar, levels, deflevel)
%----------------------------------------------------------------------
% level and lvlind will not be equal for all models, (and probably
% shouldn't for TIEGCM ... but for future expansion ...
if (nargin == 3), lvlind = deflevel; else lvlind = 1; end

fprintf('Default level (index) is  %d, if this is OK, <cr>;\n',lvlind)
fprintf('If not, enter a level between %d and %d, inclusive ...\n', ...
                         1,length(levels))
varstring = input('we''ll use the closest (no syntax required)\n','s');

if ~isempty(varstring), lvlind = str2num(varstring); end

level  = levels(lvlind);



function [lon, lonind] = GetLongitude(pgvar, lons, deflon)
%----------------------------------------------------------------------
if (nargin == 3), lon = deflon; else lon = 255.0-360.0; end

fprintf('Default longitude is %f, if this is OK, <cr>;\n',lon)
fprintf('If not, enter a longitude between %.2f and %.2f, we use the closest.\n', ...
                         min(lons),max(lons))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lon  = str2num(varstring); end

d      = abs(lon - lons);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible
lonind = ind(1);             % use the first one
lon    = lons(lonind);



function [lat, latind] = GetLatitude(pgvar, lats, deflat)
%----------------------------------------------------------------------
if (nargin == 3), lat = deflat; else lat = 40.0; end

fprintf('Default latitude is %f, if this is OK, <cr>;\n',lat)
fprintf('If not, enter a latitude between %.2f and %.2f, we use the closest.\n', ...
                         min(lats),max(lats))
varstring = input('(no syntax required)\n','s');

if ~isempty(varstring), lat = str2num(varstring); end

d      = abs(lat - lats);    % crude distance
ind    = find(min(d) == d);  % multiple minima possible
latind = ind(1);             % use the first one
lat    = lats(latind);



function varexist(filename, varnames)
%% We already know the file exists by this point.
% Lets check to make sure that file contains all needed variables.
% Fatally die if they do not exist.

nvars  = length(varnames);
gotone = ones(1,nvars);

for i = 1:nvars
   gotone(i) = nc_isvar(filename,varnames{i});
   if ( ~ gotone(i) )
      fprintf('\n%s is not a variable in %s\n',varnames{i},filename)
   end
end

if ~ all(gotone)
   error('missing required variable ... exiting')
end

