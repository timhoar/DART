function vars = CheckModel(fname)
%% CheckModel   tries to ensure that a netcdf file has what we expect. 
%
% vars is a structure containing a minimal amount of metadata about the netCDF file.
% 
% EXAMPLE:
% fname = 'Prior_Diag.nc';
% vars = CheckModel(fname) 

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if ( exist(fname,'file') ~= 2 ), error('%s does not exist.',fname); end

% Get some information from the file
model  = nc_attget(fname,nc_global,'model');
num_copies = dim_length(fname,'copy'); % determine # of ensemble members
[ens_size, ens_indices] = get_ensemble_indices(fname);
times      = nc_varget(fname,'time');
timeunits  = nc_attget(fname,'time','units');
timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin = datenum(timebase(1),timebase(2),timebase(3));
dates      = times + timeorigin;
num_times  = length(dates);

clear times timeunits timebase timeorigin

if (isempty(model)) 
   error('%s has no ''model'' global attribute.',fname)
end

switch lower(model)

   case {'9var','lorenz_63','lorenz_84','ikeda'}

      num_vars      = dim_length(fname,'StateVariable'); % determine # of state varbls
      StateVariable =  nc_varget(fname,'StateVariable');

      def_state_vars = zeros(1,num_vars);    % for use as a subscript array, 
      def_state_vars(:) = StateVariable(:);  % def_state_vars must be a row vector.

      vars = struct('model',model, ...
              'def_var','state', ...
              'num_state_vars',num_vars, ...
              'num_copies',num_copies, ...
              'num_ens_members',ens_size, ...
              'ensemble_indices',ens_indices, ...
              'time',dates, ...
              'time_series_length',num_times, ...
              'min_state_var',min(StateVariable), ...
              'max_state_var',max(StateVariable), ...
              'def_state_vars',def_state_vars);

      vars.fname = fname;

   case {'lorenz_96', 'lorenz_04'}

      num_vars      = dim_length(fname,'StateVariable'); % determine # of state varbls
      StateVariable =  nc_varget(fname,'StateVariable');

      % The only trick is to pick an equally-spaced subset of state 
      % variables for the default.

      def_state_vars = round([1 , num_vars/3 , 2*num_vars/3]);

      vars = struct('model',model, ...
              'def_var','state', ...
              'num_state_vars',num_vars, ...
              'num_copies',num_copies, ...
              'num_ens_members',ens_size, ...
              'ensemble_indices',ens_indices, ...
              'time',dates, ...
              'time_series_length',num_times, ...
              'min_state_var',min(StateVariable), ...
              'max_state_var',max(StateVariable), ...
              'def_state_vars',def_state_vars);

      vars.fname = fname;

   case 'forced_lorenz_96'

      % This model has the state variables replicated, so there is a difference
      % between num_state_vars and the length of the state variable.
      forcing           = nc_attget(fname, nc_global, 'model_forcing');
      delta_t           = nc_attget(fname, nc_global, 'model_delta_t');
      time_step_days    = nc_attget(fname, nc_global, 'model_time_step_days');
      time_step_seconds = nc_attget(fname, nc_global, 'model_time_step_seconds');
      num_model_vars    = nc_attget(fname, nc_global, 'model_num_state_vars');

      num_vars      = dim_length(fname,'StateVariable'); % determine # of state varbls

      % The only trick is to pick an equally-spaced subset of state 
      % variables for the default.

      def_state_vars = round([1 , num_model_vars/3 , 2*num_model_vars/3]);
      def_force_vars = num_model_vars + def_state_vars;

      vars = struct('model',model, ...
              'def_var','state', ...
              'num_state_vars',num_vars, ...
              'num_model_vars',num_model_vars, ...
              'num_force_vars',num_vars - num_model_vars, ...
              'num_copies',num_copies, ...
              'num_ens_members',ens_size, ...
              'ensemble_indices',ens_indices, ...
              'time',dates, ...
              'time_series_length',num_times, ...
              'min_state_var',      1, 'max_state_var', num_vars, ...
              'min_model_var',      1, 'max_model_var', num_model_vars, ...
              'min_force_var',      1, 'max_force_var', num_vars - num_model_vars, ...
              'def_state_vars',def_state_vars, ...
              'def_force_vars',def_force_vars, ...
              'forcing',forcing, ...
              'delta_t',delta_t, ...
              'time_step_days', time_step_days, ...
              'time_step_seconds', time_step_seconds);

      vars.fname = fname;

   case 'lorenz_96_2scale'

      num_X  = dim_length(fname,'Xdim'); % # of X variables
      Xdim   =  nc_varget(fname,'Xdim');

      num_Y  = dim_length(fname,'Ydim'); % # of Y variables
      Ydim   =  nc_varget(fname,'Ydim');

      % The only trick is to pick an equally-spaced subset of state 
      % variables for the default.

      def_X_inds = round([1 , num_X/3 , 2*num_X/3]);
      def_Y_inds = round([1 , num_Y/3 , 2*num_Y/3]);

      vars = struct('model',model, ...
              'def_var','X', ...
              'num_state_vars',num_X, ...
              'num_copies',num_copies, ...
              'num_ens_members',ens_size, ...
              'ensemble_indices',ens_indices, ...
              'time',dates, ...
              'time_series_length',num_times, ...
              'min_state_var',min(Xdim), 'max_state_var',max(Xdim), ...
              'min_X_var',    min(Xdim), 'max_X_var',    max(Xdim), ...
              'min_Y_var',    min(Ydim), 'max_Y_var',    max(Ydim), ...
              'def_state_vars',def_X_inds, ...
              'def_Y_inds', def_Y_inds);

      vars.fname = fname;

   case 'simple_advection'

      num_locs = dim_length(fname,'loc1d'); % # of X variables
      loc1d    =  nc_varget(fname,'loc1d');

      if ( nc_isvar(fname,'state') )
         varnames = {'state'};
         def_inds = [1 13 27];
      else
         varnames = get_DARTvars(fname);
         def_inds = round([1 , num_locs/3 , 2*num_locs/3]);
      end

      vars = struct('model'       ,model, ...
              'loc1d'             ,loc1d, ...
              'num_copies'        ,num_copies, ...
              'num_ens_members'   ,ens_size, ...
              'ensemble_indices'  ,ens_indices, ...
              'time'              ,dates, ...
              'time_series_length',num_times, ...
              'model_size'        ,length(varnames)*length(loc1d), ...
              'def_var'           ,varnames{1}, ...
              'min_state_var'     ,1, ...
              'max_state_var'     ,num_locs, ...
              'num_state_vars'    ,num_locs, ...
              'def_state_vars'    ,def_inds, ...
              'num_vars'          ,length(varnames));

      vars.vars  = varnames;
      vars.fname = fname;

   case 'wrf'

      % requires a 'domain' and 'bottom_top_d01' dimension.
      % without both of these, it will fail in an ugly fashion.

      varnames = get_DARTvars(fname);
      num_vars = length(varnames);
      num_domains = dim_length(fname,'domain');
      num_levels  = dim_length(fname,'bottom_top_d01');

      vars = struct('model',model, ...
              'def_var',varnames{1}, ...
              'def_state_vars',[], ...
              'num_state_vars',num_vars, ...
              'num_copies',num_copies, ...
              'num_ens_members',ens_size, ...
              'ensemble_indices',ens_indices, ...
              'time',dates, ...
              'time_series_length',num_times, ...
              'num_unstaggered_levels',num_levels, ...
              'num_domains',num_domains);

      vars.vars  = varnames;
      vars.fname = fname;

   case {'clm'}

      varnames = get_DARTvars(fname);
      num_vars = length(varnames);
      dinfo    = nc_getdiminfo(fname,'levtot');
      nlevels  = dinfo.Length;

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy) );

      vars.vars  = varnames;
      vars.fname = fname;

   case {'cam','tiegcm','fms_bgrid','pe2lyr','mitgcm_ocean','pbl_1d','mpas_atm'}

      varnames    = get_DARTvars(fname);
      num_vars    = length(varnames);

      vars = struct('model',model, ...
              'def_var',varnames{1}, ...
              'def_state_vars',[], ...
              'num_state_vars',num_vars, ...
              'num_copies',num_copies, ...
              'num_ens_members',ens_size, ...
              'ensemble_indices',ens_indices, ...
              'time',dates, ...
              'time_series_length',num_times);

      vars.vars  = varnames;
      vars.fname = fname;

   case 'sqg'

      % Since this is a 3D model, only the most rudimentary information
      % is gathered here. Each plot requires different information,
      % so there is a separate function (GetSqgInfo.m) that gets
      % the information for each specific plot type.

      varnames = {'theta'};
      num_vars = length(varnames);

      vars = struct('model',model, ...
              'num_state_vars',num_vars, ...
              'num_ens_members',num_copies, ...
              'time_series_length',num_times, ...
              'min_ens_mem',min(copy), ...
              'max_ens_mem',max(copy) );

      vars.vars  = varnames;
      vars.fname = fname;

   otherwise

      error('model %s unknown',model)

end


function x = dim_length(fname,dimname)
% Check for the existence of the named dimension and return it
% if it exists. If it does not, error out with a useful message.

info = nc_info(fname);
n    = length(dimname);
x    = [];
for i = 1:length(info.Dimension),
   if ( strncmp(info.Dimension(i).Name, dimname, n) > 0 )
      x = info.Dimension(i).Length;
      break
   end
end

if isempty(x)
   error('%s has no dimension named %s',fname,dimname)
end
