function PlotPhaseSpace( pinfo )
%% PlotPhaseSpace: Plots trajectories of 3 variables for any ensemble member
%
%
% PlotPhaseSpace is intended to be called by 'plot_phase_space'
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: PlotPhaseSpace( pinfo )
%
% STRUCTURE COMPONENTS FOR low-order models
% fname        name of netCDF DART file
% ens_mem      ensemble member metadata string
% var1name     name of state variable to plot as 'X'
% var1ind      ID   of state variable to plot as 'X'
% var2name     name of state variable to plot as 'Y'
% var2ind      ID   of state variable to plot as 'Y'
% var3name     name of state variable to plot as 'Z'
% var3ind      ID   of state variable to plot as 'Z'
% ltype        line type (see 'help plot' for details)
%
% Example 1   ( 9 variable model )
%%--------------------------------------------------------
% pinfo.fname      = 'True_State.nc';
% pinfo.ens_mem    = 'true state';    % true state only has 1 ens mem ...
% pinfo.var1name   = 'state';         % 9var netCDF has only 1 flavor variable
% pinfo.var2name   = 'state';
% pinfo.var3name   = 'state';
% pinfo.var1ind    = 3;
% pinfo.var2ind    = 6;
% pinfo.var3ind    = 7;
% pinfo.ltype      = 'b-'; % solid blue line
% PlotPhaseSpace( pinfo )
%
% that worked so well, lets overlay another (using the same state variables)
%
% hold on;
% pinfo.fname      = 'Prior_Diag.nc';
% pinfo.ens_mem    = 'ensemble member4';               % why not?
% pinfo.ltype      = 'r:';            % plot it in a red 'dotted' line
% PlotPhaseSpace( pinfo )
%
% note the legend has both lines annotated.

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if ( exist(pinfo.fname,'file') ~= 2 ), error('file %s does not exist.',pinfo.fname), end

model = nc_attget(pinfo.fname, nc_global, 'model');

% Get some information for the 'X','Y' variable
[X.num_times, X.num_copies, X.num_vars] = parse_varshape(pinfo.fname, pinfo.var1name);
[Y.num_times, Y.num_copies, Y.num_vars] = parse_varshape(pinfo.fname, pinfo.var2name);

if ( isfield( pinfo,'var3name') ) % Get some information for the 'Z' variable
[Z.num_times, Z.num_copies, Z.num_vars] = parse_varshape(pinfo.fname, pinfo.var3name);
end

switch lower(model)

   case {'ikeda'}   % only two state variables

      ens_mem_id = get_copy_index(pinfo.fname, pinfo.ens_mem);  % errors out if no ens_mem

      x = get_var_series(pinfo.fname, pinfo.var1name, ens_mem_id, pinfo.var1ind);
      y = get_var_series(pinfo.fname, pinfo.var2name, ens_mem_id, pinfo.var2ind);

      h = plot(x,y,pinfo.ltype);
      axis image

      % datmin = min(min(x,y));
      % datmax = max(max(x,y));
      % axis([datmin datmax datmin datmax])

      % If there is no legend, we are assuming this is the first plot on this
      % axis. We need to add a title, legend, axes labels, etc.
      [legh, objh, outh, outm] = legend;
      if (isempty(legh))

         % title(sprintf('%s ensemble member %d of %s',model,ens_mem_id,pinfo.fname), ...
         %    'interpreter','none','fontweight','bold')
         title('The Ikeda phase space','fontweight','bold')
         xlabel(sprintf('%s variable # %d',pinfo.var1name, pinfo.var1ind))
         ylabel(sprintf('%s variable # %d',pinfo.var2name, pinfo.var2ind))

         s = sprintf('%d %d %s %s %s', pinfo.var1ind, pinfo.var2ind, ...
                     pinfo.fname, model, pinfo.ens_mem);
         h = legend(s,0);
         [legh, objh, outh, outm] = legend;
         set(objh(1),'interpreter','none')

      else
         % Must add salient information to the legend.
         % legh     handle to the legend axes
         % objh     handle for the text, lines, and patches in the legend
         % outh     handle for the lines and patches in the plot
         % outm     cell array for the text in the legend
         nlines = length(outm);
         outm{nlines+1} = sprintf('%d %d %s %s %s', ...
              pinfo.var1ind, pinfo.var2ind, pinfo.fname, ...
              model, pinfo.ens_mem);
         [legh, objh, outh, outm] = legend([outh; h],outm,0);

         set(objh(1:nlines+1),'interpreter','none')
      end
      legend boxoff

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','simple_advection'}

      BulletProof(pinfo,X,Y,Z)          % rudimentary bulletproofing

      ens_mem_id = get_copy_index(pinfo.fname, pinfo.ens_mem);  % errors out if no ens_mem

      x = get_var_series(pinfo.fname, pinfo.var1name, ens_mem_id, pinfo.var1ind);
      y = get_var_series(pinfo.fname, pinfo.var2name, ens_mem_id, pinfo.var2ind);
      z = get_var_series(pinfo.fname, pinfo.var3name, ens_mem_id, pinfo.var3ind);

      % There is no model-dependent segment ...
      % As long as you have three variables, this works for all models.

      h = plot3(x,y,z,pinfo.ltype);

      % If there is no legend, we are assuming this is the first plot on this
      % axis. We need to add a title, legend, axes labels, etc.
      [legh, objh, outh, outm] = legend;
      if (isempty(legh))

         % title(sprintf('%s ensemble member %d of %s',model,ens_mem_id,pinfo.fname), ...
         %    'interpreter','none','fontweight','bold')
         title('The Attractor','fontweight','bold')
         xlabel(sprintf('%s variable # %d',pinfo.var1name, pinfo.var1ind))
         ylabel(sprintf('%s variable # %d',pinfo.var2name, pinfo.var2ind))
         zlabel(sprintf('%s variable # %d',pinfo.var3name, pinfo.var3ind))

         s = sprintf('%d %d %d %s %s %s', pinfo.var1ind, pinfo.var2ind, ...
                     pinfo.var3ind, pinfo.fname, model, pinfo.ens_mem);
         h = legend(s,0);
         [legh, objh, outh, outm] = legend;
         set(objh(1),'interpreter','none')

      else
         % Must add salient information to the legend.
         % legh     handle to the legend axes
         % objh     handle for the text, lines, and patches in the legend
         % outh     handle for the lines and patches in the plot
         % outm     cell array for the text in the legend
         nlines = length(outm);
         outm{nlines+1} = sprintf('%d %d %d %s %s %s', ...
              pinfo.var1ind, pinfo.var2ind, pinfo.var3ind, pinfo.fname, ...
              model, pinfo.ens_mem);
         [legh, objh, outh, outm] = legend([outh; h],outm,0);

         set(objh(1:nlines+1),'interpreter','none')
      end
      legend boxoff

   case {'fms_bgrid','pe2lyr','mitgcm_ocean','wrf','cam'}

      ens_mem_id = get_copy_index(pinfo.fname, pinfo.ens_mem);   % errors out if no ens_mem

      x = Get1Copy(pinfo.fname, pinfo.var1name, ens_mem_id,  ...
                  pinfo.var1_lvlind, pinfo.var1_latind, pinfo.var1_lonind);
      y = Get1Copy(pinfo.fname, pinfo.var2name, ens_mem_id,  ...
                  pinfo.var2_lvlind, pinfo.var2_latind, pinfo.var2_lonind);
      z = Get1Copy(pinfo.fname, pinfo.var3name, ens_mem_id,  ...
                  pinfo.var3_lvlind, pinfo.var3_latind, pinfo.var3_lonind);

      % There is no model-dependent segment ...
      % As long as you have three variables, this works for all models.

      h = plot3(x,y,z,pinfo.ltype);

      % If there is no legend, we are assuming this is the first plot on this
      % axis. We need to add a title, legend, axes labels, etc.
      [legh, objh, outh, outm] = legend;
      if (isempty(legh))

         % title(sprintf('%s ensemble member %d of %s',model,ens_mem_id,pinfo.fname), ...
         %    'interpreter','none','fontweight','bold')
         title('The Attractor','fontweight','bold')
         hx = xlabel(sprintf('%s %d %.2f %.2f', ...
               pinfo.var1name, pinfo.var1_lvl, pinfo.var1_lat, pinfo.var1_lon));
         hy = ylabel(sprintf('%s %d %.2f %.2f', ...
               pinfo.var2name, pinfo.var2_lvl, pinfo.var2_lat, pinfo.var2_lon));
         hz = zlabel(sprintf('%s %d %.2f %.2f', ...
               pinfo.var3name, pinfo.var3_lvl, pinfo.var3_lat, pinfo.var3_lon));
         set(hx,'interpreter','none')
         set(hy,'interpreter','none')
         set(hz,'interpreter','none')

         s = sprintf('%s %s %s %s %s %s', pinfo.var1name, pinfo.var2name, pinfo.var3name, ...
                                              model, pinfo.fname, pinfo.ens_mem);
         h = legend(s,0);
         [legh, objh, outh, outm] = legend;
         set(objh(1),'interpreter','none')

      else
         % Must add salient information to the legend.
         % legh     handle to the legend axes
         % objh     handle for the text, lines, and patches in the legend
         % outh     handle for the lines and patches in the plot
         % outm     cell array for the text in the legend
         nlines = length(outm);
         outm{nlines+1} = sprintf('%s %s %s %s %s %s', pinfo.var1name, ...
               pinfo.var2name, pinfo.var3name, model, pinfo.fname, pinfo.ens_mem);
         [legh, objh, outh, outm] = legend([outh; h],outm,0);

         set(objh(1:nlines+1),'interpreter','none')
      end
      legend boxoff

   otherwise

      error('model %s not implemented yet', model)

end


%======================================================================
% Subfunctions
%======================================================================



function [nT, nC, n3] = parse_varshape(fname,varname)

y = nc_isvar(fname,varname);
if (y < 1)
   error('%s has no variable named %s ',fname,varname)
end

nt = 0;
nC = 0;
n3 = 0;

varinfo = nc_getvarinfo(fname,varname);

% for i = 1:length(varinfo.Dimension)
for i = 1:3  % only want/need the first 3 dimensions.
   switch( lower(varinfo.Dimension{i}))
      case 'time'
         nT = varinfo.Size(i);
      case 'copy'
         nC = varinfo.Size(i);
      case 'StateVariable'
         n3 = varinfo.Size(i);
      otherwise
         n3 = varinfo.Size(i);
   end
end


function var = Get1Copy(fname, varname, copyindex, lvlind, latind, lonind)

% Gets a time-series of a single specified copy of a prognostic variable
% at a particular 3D location (level, lat, lon)

myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
myinfo.levelindex = lvlind;
myinfo.latindex   = latind;
myinfo.lonindex   = lonind;
[start, count]    = GetNCindices(myinfo,'diagn',varname);

var = nc_varget(fname, varname, start, count); % 'bob' is only 2D




function BulletProof(pinfo,X,Y,Z)

if ( (pinfo.var1ind > X.num_vars) | (pinfo.var1ind < 1) )
   fprintf('\n%s has %d %s variables\n',pinfo.fname,X.num_vars,pinfo.var1name)
   fprintf('%d  <= ''var1'' <= %d\n',1,X.num_vars)
   error('var1 (%d) out of range',pinfo.var1ind)
end

if ( (pinfo.var2ind > Y.num_vars) | (pinfo.var2ind < 1) )
   fprintf('\n%s has %d %s variables\n',pinfo.fname,Y.num_vars,pinfo.var2name)
   fprintf('%d  <= ''var2'' <= %d\n',1,Y.num_vars)
   error('var2 (%d) out of range',pinfo.var2ind)
end

if ( (pinfo.var3ind > Z.num_vars) | (pinfo.var3ind < 1) )
   fprintf('\n%s has %d %s variables\n',pinfo.fname,Z.num_vars,pinfo.var3name)
   fprintf('%d  <= ''var3'' <= %d\n',1,Z.num_vars)
   error('var3 (%d) out of range',pinfo.var3ind)
end
