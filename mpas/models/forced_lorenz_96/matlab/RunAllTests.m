function RunAllTests(dummy)
%% RunAllTests.m

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (nargin() > 0)
   interactive = 1;
else
   interactive = 0;
end

figure(1)
if (interactive)
 fprintf('Starting %s\n','plot_bins');
 clear truth_file diagn_file; close all; plot_bins
 fprintf('Finished %s ... pausing, hit any key\n','plot_bins'); pause
 fprintf('Starting %s\n','plot_ens_err_spread');
 plot_ens_err_spread
 fprintf('Finished %s ... pausing, hit any key\n','plot_ens_err_spread'); pause
 fprintf('Starting %s\n','plot_ens_time_series');
 plot_ens_time_series
 fprintf('Finished %s ... pausing, hit any key\n','plot_ens_time_series'); pause
 fprintf('Starting %s\n','plot_ens_mean_time_series');
 plot_ens_mean_time_series
 fprintf('Finished %s ... pausing, hit any key\n','plot_ens_mean_time_series'); pause
end

 fprintf('Starting %s\n','PlotBins');
 clear pinfo; close all;

 pinfo = CheckModelCompatibility('True_State.nc','Prior_Diag.nc');
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.diagn_file);
 pinfo.var            = 'state';
 pinfo.var_inds       = [1 13 27];

 PlotBins(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotBins'); pause

 fprintf('Starting %s\n','PlotEnsErrSpread');
 close all; PlotEnsErrSpread(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotEnsErrSpread'); pause

 fprintf('Starting %s\n','PlotEnsTimeSeries');
 close all; PlotEnsTimeSeries(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotEnsTimeSeries'); pause

 fprintf('Starting %s\n','PlotEnsMeanTimeSeries');
 close all; PlotEnsMeanTimeSeries(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotEnsMeanTimeSeries'); pause

%% ----------------------------------------------------------
% plot_correl
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n', 'plot_correl');
 clear diagn_file; close all; plot_correl
 fprintf('Finished %s ... pausing, hit any key\n','plot_correl'); pause
end

 fprintf('Starting %s\n','PlotCorrel');
 clear pinfo; clf

 pinfo                    = CheckModel('Prior_Diag.nc');
 pinfo.base_var           = 'state';
 pinfo.base_var_index     =  30;
 pinfo.base_time          =  300;

 PlotCorrel(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotCorrel'); pause

%% -----------------------------------------------------------
% plot_phase_space
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_phase_space');
 clear fname; close all; plot_phase_space
 fprintf('Finished %s ... pausing, hit any key\n','plot_phase_space'); pause
end

 fprintf('Starting %s\n','PlotPhaseSpace');
 clear pinfo; clf

 pinfo.fname    = 'True_State.nc';
 pinfo.model    = 'Forced_Lorenz_96';
 pinfo.var1name = 'state';
 pinfo.var2name = 'state';
 pinfo.var3name = 'state';
 pinfo.var1ind  = 10;
 pinfo.var2ind  = 20;
 pinfo.var3ind  = 30;
 pinfo.ens_mem  = 'true state';
 pinfo.ltype    = 'k-';

 PlotPhaseSpace(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotPhaseSpace'); pause

%% ----------------------------------------------------------
% plot_reg_factor
%------------------------------------------------------------
% plot_reg_factor

%% ----------------------------------------------------------
% plot_sawtooth
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_sawtooth');
 clear truth_file posterior_file prior_file; close all; plot_sawtooth
 fprintf('Finished %s ... pausing, hit any key\n','plot_sawtooth'); pause
end

 fprintf('Starting %s\n','PlotSawtooth');
 clear pinfo; close all

 pinfo    = CheckModelCompatibility('Prior_Diag.nc','Posterior_Diag.nc');
 pinfo.prior_time     = pinfo.truth_time;
 pinfo.prior_file     = pinfo.truth_file;
 pinfo.posterior_time = pinfo.diagn_time;
 pinfo.posterior_file = pinfo.diagn_file;
 pinfo.truth_file     = 'True_State.nc';
 pinfo = rmfield(pinfo,{'diagn_file','truth_time','diagn_time'});
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.prior_file);
 pinfo.var            = 'state';
 pinfo.var_inds       = [1 20 40];
 pinfo.copyindices    = [7 12 17];

 PlotSawtooth(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotSawtooth'); pause

%% -----------------------------------------------------------
% plot_smoother_err
%------------------------------------------------------------
% plot_smoother_err

%% ----------------------------------------------------------
% plot_total_err
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_total_err');
 clear truth_file diagn_file; close all; plot_total_err
 fprintf('Finished %s ... pausing, hit any key\n','plot_total_err'); pause
end

 fprintf('Starting %s\n','PlotTotalErr');
 clear pinfo; clf

 pinfo = CheckModelCompatibility('True_State.nc','Prior_Diag.nc');

 PlotTotalErr(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotTotalErr'); pause

%% ----------------------------------------------------------
% plot_var_var_correl
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_var_var_correl');
 clear fname; close all; plot_var_var_correl
 fprintf('Finished %s ... pausing, hit any key\n','plot_var_var_correl'); pause
end

 fprintf('Starting %s\n','PlotVarVarCorrel');
 clear pinfo; clf

 pinfo                 = CheckModel('Prior_Diag.nc');
 pinfo.base_var        = pinfo.def_var;
 pinfo.state_var       = pinfo.def_var;
 pinfo.base_var_index  = 30;
 pinfo.base_time       = 300;
 pinfo.state_var_index = 60;

 PlotVarVarCorrel(pinfo)
 fprintf('Finished %s ... pausing, hit any key\n','PlotVarVarCorrel'); pause

%% ----------------------------------------------------------
% plot_jeff_correl - correlation evolution
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_jeff_correl');
 clear fname; close all; plot_jeff_correl
 fprintf('Finished %s ... pausing, hit any key\n','plot_jeff_correl'); pause
end

 fprintf('Starting %s\n','PlotJeffCorrel');
 PlotJeffCorrel(pinfo)
 fprintf('Finished %s\n','PlotJeffCorrel')
