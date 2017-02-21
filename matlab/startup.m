%% startup.m  IFF a $HOME/matlab/startup.m exists, it is executed automatically at matlab's startup.
%
% The netcdf toolbox is needed for any/all DART matlab diagnostics, so this
% block tries to locate that particular startup script.
% The beauty of addpath is that if the desired directory is already
% in your path, nothing happens, so there is no harm trying.

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%% Try to intelligently add the general DART tools.

mydir      = pwd;
dartloc    = strfind(mydir,'/models/'  )-1;
dartlabloc = strfind(mydir,'/DART_LAB/')-1;
if (isempty(dartloc) && isempty(dartlabloc))
   return
end

fprintf('\nWelcome to DART ...\n')
fprintf('\nYour current directory is  %s\n',mydir)

% DART/matlab directory ...

dartpath = sprintf('%s/matlab',mydir(1:dartloc));
if (exist(dartpath,'dir')==7)
   addpath(dartpath,'-BEGIN');
   fprintf('Using general tools in     %s\n',dartpath)
end

% add the observation-space DART tools ...

dartpath = sprintf('%s/diagnostics/matlab',mydir(1:dartloc));
if ( ~isempty(dartloc) )
   addpath(dartpath,'-BEGIN');
   fprintf('observation-space tools in %s\n',dartpath)
end

% add the DART_LAB/matlab directory ...

dartpath = sprintf('%s/DART_LAB/matlab',mydir(1:dartloc));
if (exist(dartpath,'dir')==7)
   addpath(dartpath,'-BEGIN');
   fprintf('Using DART_LAB tools in    %s\n',dartpath)
end

% Try to intelligently add the DART model-specific tools.
% If the cwd is a '<model>/work' directory, check to see if there is a
% parallel '<model>/matlab' directory.

dartloc  = strfind(mydir,'/work')-1;
if ( ~isempty(dartloc) )
   dartpath = sprintf('%s/matlab',mydir(1:dartloc));
   if (exist(dartpath,'dir') == 7 )
      addpath(dartpath,'-BEGIN');
      fprintf('model-specific scripts in  %s\n',dartpath)
   end
end

% summarize

truth_file = fullfile(mydir,'perfect_output.nc');
diagn_file = fullfile(mydir,'preassim.nc');

disp(' ')
fprintf('the default data directory is          %s\n',mydir)
fprintf('which means your default TRUTH file is %s\n',truth_file)
fprintf('and your default    DIAGNOSTIC file is %s\n',diagn_file)
disp('To change your defaults, set ''truth_file'' and/or ''diagn_file'' accordingly.')

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
