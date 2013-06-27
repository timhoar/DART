%% DART:script4  With ensemble methods, prior is available only as a 'random' sample

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Begin by setting random seed to a nice controlled
% initial value that gives nice plots
randn('state', 0);

% Generate a random sample of the prior
x_prior = normrnd(-1.0, 1.2, 5, 1);
% Adjust this to have an exact mean and variance consistent with desired distribution
x_prior = (x_prior - mean(x_prior)) * 1.2 / std(x_prior) - 1.0

% Plot the prior distribution
y_prior = [0.05, 0.05, 0.05, 0.05, 0.05];
h_plot = plot(x_prior, y_prior, 'g*')
set(h_plot, 'markersize', 18)
set(h_plot, 'linewidth', 3)
set(h_plot, 'color', [0 0.73 0]);

axis([-4 4 0 0.6]);
set(gca, 'fontsize', 24);
text(-2.0, 0.1, 'Prior Ensemble', 'fontsize', 24);
% Set the size and position of the figure box on the screen
set(gcf, 'Units', 'Inches');
set(gcf, 'position', [1 1 10.5 5.0]);
set(gca, 'linewidth', 2);
grid on;
% Set the shape of the plot box
pbaspect([2.4 1 1]);
ylabel('Probability');
hold on;
pause

% Setup the printing characteristics
% PRINTING IS SCREWED UP WITH CURRENT MATLAB VERSION, DO MANUALLY                        
%set(gca, 'xtickmode', 'manual', 'ytickmode', 'manual');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
print -depsc s04f01.eps;                       

% Plot a gaussian fit to ensemble
x = -5:0.01:5;
prior = normpdf(x, -1.0, 1.2);
obs = normpdf(x, 1.0, 0.8);
product = prior .* obs;
hhh = plot(x, prior, 'g', 'linewidth', 3)
set(hhh, 'color', [0 0.73 0]);
text(-3.0, 0.3, 'Prior PDF', 'fontsize', 24);
pause;
print -depsc s04f02.eps;                       

plot(x, obs, 'r', 'linewidth', 3)
text(1.6, 0.4, 'Obs. Likelihood', 'fontsize', 24);
pause;
print -depsc s04f03.eps;                       

% Need to get integrated value of the observed value, y2
pb = (sum(product) * 0.01); 
posterior = product./pb;
plot(x, posterior, 'b', 'linewidth', 3)
text(-1.9, 0.55, 'Posterior PDF', 'fontsize', 24);
print -depsc s04f04.eps;                       

