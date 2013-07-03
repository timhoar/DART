% startup.m  IFF a $HOME/matlab/startup.m exists, it is executed automatically at matlab's startup.


% The netcdf toolbox is needed for any/all DART matlab diagnostics, so this
% block tries to locate that particular startup script.
% The beauty of addpath is that if the desired directory is already 
% in your path, nothing happens, so there is no harm trying.

if ( exist('/contrib/matlab/ncstartup.m') == 2 ) 
   addpath  /contrib/matlab
   ncstartup;                     % Adds the netCDF operators
elseif ( exist('/usr/local/matlab/ncstartup.m') == 2 ) 
   addpath      /usr/local/matlab
   ncstartup;                     % Adds the netCDF operators
end

% See if we have succeeded in adding the netcdf operators.

if ( exist('getnc') ~= 2 ) 
   disp('Sorry. Unable to locate the netCDF matlab operators.')
   error('The DART diagnostics will not run.')
end

% Try to intelligently add the general DART tools.

mydir    = pwd;
dartloc  = strfind(mydir,'/DART/')+4;
dartpath = sprintf('%s/matlab',mydir(1:dartloc));

disp(sprintf('\nWelcome to DART ...'))
disp(sprintf('\nYour current directory is  %s',mydir))

if ( ~isempty(dartloc) ) 
   path(dartpath,path);
   disp(sprintf('Using general tools in     %s',dartpath))
end

% Try to intelligently add the observation-space DART tools.

dartpath = sprintf('%s/diagnostics/matlab',mydir(1:dartloc));

if ( ~isempty(dartloc) ) 
   path(dartpath,path);
   disp(sprintf('observation-space tools in %s',dartpath))
end

% Try to intelligently add the DART model-specific tools.
% If the cwd is a '<model>/work' directory, check to see if there is a 
% parallel '<model>/matlab' directory.

mydir    = pwd;
dartloc  = strfind(mydir,'/work')-1;
dartpath = sprintf('%s/matlab',mydir(1:dartloc));

if ( ~isempty(dartloc) ) 
   if ( exist(dartpath,'dir') == 7 )
      path(dartpath,path);
      disp(sprintf('Using matlab scripts in     %s',dartpath))
   end
end

% Customizations specific for DART:

datadir = '.';
truth_file = fullfile(datadir,'True_State.nc');
diagn_file = fullfile(datadir,'Prior_Diag.nc');

disp(' ')
disp(sprintf('the default data directory is          %s',datadir))
disp(sprintf('which means your default TRUTH file is %s',truth_file))
disp(sprintf('and your default    DIAGNOSTIC file is %s',diagn_file))
disp('To change your defaults, set ''truth_file'' and/or ''diagn_file'' accordingly.')
