function check_forecast(component,fname)
%% check_forecast is meant to assist verifying the output of obs_seq_coverage.f90 and obs_seq_to_stations.f90
%
% To check the output of obs_seq_coverage.f90:
% fname     = 'obsdef_mask.nc';
% check_forecast('coverage',fname)
%
% To check the output of obs_seq_to_stations.f90:
% fname     = 'METAR_U_10_METER_WIND_forecast.nc';
% check_forecast('forecast',fname)
%

%% DART software - Copyright © 2004 - 2010 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

switch lower(component)
   case {'forecast','stations'}
      check_forecast_component(fname);
   otherwise
      check_coverage_component(fname);
end

function check_forecast_component(fname)

indx        = strfind(fname,'_forecast.nc');
varname     = fname(1:indx-1);

X           = nc_varget(fname,varname);
analysisT   = nc_varget(fname,'analysisT');
copy        = nc_varget(fname,'copy');
stations    = nc_varget(fname,'stations');
levels      = nc_varget(fname,'levels');
nmembers    = nc_varget(fname,'nmembers');
fcstlead    = nc_varget(fname,'forecast_lead');
CopyStrings = nc_varget(fname,'CopyMetaData');
qc          = nc_varget(fname,'original_qc');
dartqc      = nc_varget(fname,'dart_qc');

xax         = [0 max(nmembers)];

fprintf('There are %d stations, %d of which are not missing the first timestep.\n', ...
         length(stations), sum(isfinite(X(:,1,1,1))) )

% X has shape (roughly) :
% 1000           3          50           5
% stations    copies      ens_size   forecast_leads

% copy 2 should be prior ... should not be identical ...
% see about other stations

for istation = 1:length(stations)

   for itime = 1:4
      obsval = X(istation,1,1,itime);  % same for all ens mems
      errvar = X(istation,3,1,itime);  % same for all ens mems

      figure(1)

      string1 = sprintf('loc %d %s qc %d dartqc %d', ...
          istation,deblank(CopyStrings(2,:)),qc(istation,itime),dartqc(istation,itime));
      subplot(2,2,itime);
      plot(nmembers,squeeze(X(istation,2,:,itime)),'o',xax,[obsval obsval],'-')
      title(string1)
      xlabel(sprintf('error variance %f',errvar))
      h = ylabel(varname);
      set(h,'Interpreter','none')

      % TJH IMPROVEME/FIXME observation value should have some sort of whisker
      % plot reflecting observation error variance

   end

   disp('pausing, hit any key ...')
   pause

end


function check_coverage_component(fname)

stationtime = nc_varget(fname,'StationTime');

stationtime(isnan(stationtime)) = 0.0;

vtime       = nc_varget(fname,'time');
timeunits   = nc_attget(fname,'time','units');
calendar    = nc_attget(fname,'time','calendar');
timebase    = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin  = datenum(timebase(1),timebase(2),timebase(3));
vtimes      = datestr(vtime + timeorigin);

%% histogram of how many stations have how many verification times
figure(1); clf; orient landscape
ntimes      = single(nc_varget(fname,'ntimes'));
nmin        = min(ntimes);
nmax        = max(ntimes);

hist(ntimes,nmin:1:nmax);
axlims    = axis;
axlims(1) = nmin - 0.5;
axlims(2) = nmax + 0.5;
axis(axlims)

blanks = zeros(size(vtimes,1),2);
blanks(:) = ' ';


for istation = 1:length(stations)

    fprintf('\nStation %d flagged as %d\n',istation,stations(istation))

    mytime = datestr(stationtime(istation,:) + timeorigin);

    [vtimes blanks mytime]

    pause

end

