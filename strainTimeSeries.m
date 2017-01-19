%%
% Script for calculating the strain rate time series for Afines data
% This uses an external file "strain.m" which finds the strain rate tensor, eigenvectors
% and eigenvalues for the given vector field. I'm including that file at the end, it can be found
% on my github page.

clear;
close all;
tic;

str = ['Data contained in: ',pwd];
fprintf(str)

fname = 'interpedData.mat';
myVars = {'adataVelInterped','dr'};
load(fname, myVars{:});
[~, ~, ~, nFrames] = size(adataVelInterped);
strainRateSeries = zeros(nFrames,1);

for frame = 1:size(adataVelInterped,4)
    interpedVX = adataVelInterped(:,:,1,frame);
    interpedVY = adataVelInterped(:,:,2,frame);

    [strainRateTensor, strainRateEvecs, strainRateEvals] = strain(interpedVX, interpedVY, dr);
    div = sum(sum(strainRateEvals,3),4);
    strainRateSeries(frame) = mean(mean(div));
end

plot(strainRateSeries, 'bo')

Max = max(strainRateSeries);
Min = min(strainRateSeries);
Mean = mean(strainRateSeries);
Std = std(strainRateSeries);
xlim([0,500])
ylimit = max([abs(Max), abs(Min)])*1.2;
ylim([-ylimit, ylimit])
txt = ['max = ', num2str(Max), ', min = ', num2str(Min)];
txt2 = ['Mean \pm std = ', num2str(Mean), '\pm', num2str(Std)];
text(250, ylimit*0.9, txt);
text(250, ylimit*0.8, txt2);
title('Mean of $\vec{\nabla} \cdot \vec{v}$ over whole domain', 'Interpreter', 'latex')


%% 
%
% Here we see how the mean divergence of the entire frame is always around zero. This arises due to
% a combination of periodic boundary conditions and mass conservation. The mass continuity equation is given
% by $\partial_t m + \vec{\nabla} \cdot \vec{p} = 0$, and since $\partial_t m = 0$ with periodic boundaries,
% we are left with $\vec{\nabla} \cdot \vec{p} = m_{bead} \vec{\nabla} \cdot \vec{v} = 0$, showing that the
% divergence of the velocity field must be zero for periodic boundary conditions.
%
% Below is the strain.m function
% 
% <include>strain.m</include>