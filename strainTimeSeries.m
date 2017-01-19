% Script for calculating the strain rate time series for Afines data
%

clear;
close all;
tic;

fname = 'interpedData.mat';
myVars = {'adataVelInterped'};
load(fname, myVars{:});
[~, ~, ~, nFrames] = size(adataVelInterped);
strainRateSeries = zeros(nFrames,1);

for frame = 1:size(adataVelInterped,4)
    interpedVX = adataVelInterped(:,:,1,frame);
    interpedVY = adataVelInterped(:,:,2,frame);

    [strainRateTensor, strainRateEvecs, strainRateEvals] = strain(interpedVX, interpedVY, dr);
    div = sum(sum(strainRateEvals,3),4);
    strainRateSeries(frame) = mean(mean(div));

    % pcolor(xGrid,yGrid,div);
    % shading interp;
    % axis image;
    % colorbar;
    % caxis([-1 1])
    % hold on;
    % quiver(xGrid, yGrid, interpedVX, interpedVY,'k')
    % pause(0.1)
    % hold off
end


