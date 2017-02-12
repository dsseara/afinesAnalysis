%%
% Script for calculating the strain rate time series for Afines data
% This uses an external file "strain.m" which finds the strain rate tensor, eigenvectors
% and eigenvalues for the given vector field. 

clear;
close all;
tic;

str = ['Data contained in: ',pwd];
fprintf(str);

fname = 'interpedData.mat';
myVars = {'adataVelInterped','dr', 'xRange', 'yRange'};
load(fname, myVars{:});
[~, ~, ~, nFrames] = size(adataVelInterped);
strainRateSeries = zeros(nFrames,1);
seriesStats.std = [];
seriesStats.mean = [];

halfSpaceSize = min(xRange(2)-xRange(1), yRange(2)-yRange(1));
for i = 0:floor(halfSpaceSize)-2 % Want to leave at least some space in the middle..
    interiorXRange = [xRange(1)+i, xRange(2)-i];
    interiorYRange = [yRange(1)+i, yRange(2)-i];
    figure()
    for frame = 1:size(adataVelInterped,4)
        interpedVX = adataVelInterped((i+1):(end-i), (i+1):(end-i), 1, frame);
        interpedVY = adataVelInterped((i+1):(end-i), (i+1):(end-i), 2, frame);

        [strainRateTensor, strainRateEvecs, strainRateEvals] = strain(interpedVX, interpedVY, dr);
        div = sum(sum(strainRateEvals,3),4);
        strainRateSeries(frame) = mean(mean(div));
    end
    seriesStats.std = [seriesStats.std std(strainRateSeries)];
    seriesStats.mean = [seriesStats.mean mean(strainRateSeries)];
    plot((1:10:2000)./4,strainRateSeries(1:10:end),'o')
    xlabel('time (s)')
    ylabel('$ \dot{\epsilon}(t) $', 'Interpreter', 'LaTex')
    txt = ['Mean $\pm$ std = ', num2str(seriesStats.mean(end)),'$\pm$' num2str(seriesStats.std(end))];
    text(250, 0.5*max(strainRateSeries), txt, 'Interpreter','LaTeX');
    txt2 = ['X-range = ', num2str(interiorXRange(1)), '-', num2str(interiorXRange(2))];
    txt3 = [', Y-range = ', num2str(interiorYRange(1)), '-', num2str(interiorYRange(2))];
    title([txt2, txt3])

end

errorbar(0:(numel(seriesStats.mean)-1), seriesStats.mean, seriesStats.std, 'o');
hold on;
plot([0,30],[0,0])
xlim([0,30])
legend('Mean Strain Rate, \mu \pm \sigma', 'Zero Line')

%% 
% 
% These plots show the mean strain rate over time over an increasingly small area in the center of the domain of the simulation. As we can see, the average strain rate is consistently zero, despite the changing box size. Although there always does seem to be a peak at the very beginning of the simulation.