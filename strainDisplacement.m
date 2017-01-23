%% Strain measured as average displacement
% 
%

clear; close all;
tic;
cd /Users/Danny/'Google Drive'/'Murrell Lab'/AFiNeS_Analysis/tf500_npolymer_500_a_motor_density_0.5_p_motor_density_1;
fname = 'simdata.mat';
vars = {'adata', 'timestep'};
load(fname, vars{:});
L = 50;

% Enforce periodic boundary conditions
numTimeSteps = 1;

xyDisplacement = mod(diff(adata(:,1:2, 1:numTimeSteps:end),1,3) + L/2, L) - L/2;
[numBeads, dim, numFrames] = size(xyDisplacement);
dx = squeeze(xyDisplacement(:,1,:));
dy = squeeze(xyDisplacement(:,2,:)); clear xyDisplacement;

cumDR = [zeros(size(dx,1),1), sqrt(cumsum(dx,2).^2 + cumsum(dy,2).^2)]./L;

meanCumDR = mean(cumDR);
h = zeros(2,1);
h(1) = plot(meanCumDR, 'DisplayName', '$\rho_a = 0.5$, $\rho_p = 1$');
hold on

cd /Users/Danny/'Google Drive'/'Murrell Lab'/AFiNeS_Analysis/tf500_npolymer_500_a_motor_density_0.1_p_motor_density_0.1;

fname = 'simdata.mat';
vars = {'adata', 'timestep'};
load(fname, vars{:});
L = 50;

% Enforce periodic boundary conditions
numTimeSteps = 1;

xyDisplacement = mod(diff(adata(:,1:2, 1:numTimeSteps:end),1,3) + L/2, L) - L/2;
[numBeads, dim, numFrames] = size(xyDisplacement);
dx = squeeze(xyDisplacement(:,1,:));
dy = squeeze(xyDisplacement(:,2,:)); clear xyDisplacement;

cumDR = [zeros(size(dx,1),1), sqrt(cumsum(dx,2).^2 + cumsum(dy,2).^2)]./L;

meanCumDR = mean(cumDR);

h(2) = plot(meanCumDR, 'DisplayName', '$\rho_a = 0.1$, $\rho_p = 0.1$');
ylabel('$\langle|x(t) - x(0)|\rangle/L$', 'Interpreter', 'Latex')
l = legend(h);
set(l, 'Interpreter', 'Latex')
title('Strain time series')

cd ..