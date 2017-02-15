%%
% This script gets the statistics out of velocity divergence data output running divVelocity.m on a simdata.mat data file that was generated using read_data.m on the actins.txt output of an afines simulation
%
% Created by Daniel Seara at 2017/02/15 11:05

clear; close all;

fname='interpedData.mat';
myVars = {'grid', 'binParams', 'divV', 'interpParams'};
load(fname, myVars{:});
fname='simdata.mat';
myVars={'params'};
load(fname, myVars{:});

if ~exist('divStatsPlots','dir')
    mkdir('divStatsPlots')
end

% First we calculate the mean of the total velocity divergence over the whole domain and plot its time series

meanDivV = squeeze(mean(mean(divV)));
time = (params.timestep(1)+binParams.DeltaT*0.5):binParams.DeltaT:(params.timestep(end)-binParams.DeltaT*0.5); % Get the half-time steps where velocity is actually calculated
plot(time, meanDivV)
meanTimeSeries = mean(meanDivV);
stdTimeSeries  =  std(meanDivV);
xlabel('time (s)')
ylabel('$\langle \vec{\nabla} \cdot \vec{v} \rangle (t)$', 'Interpreter', 'latex')
legend(sprintf('\\mu \\pm \\sigma = %.02e \\pm %.02e', meanTimeSeries, stdTimeSeries))
title('Average over whole domain')%, 'Interpreter','none')
prettyFig;

if isunix
    saveas(gcf, [pwd, '/divStatsPlots/' , 'totalMean', '.tif'])
elseif ispc
    saveas(gcf, [pwd, '\divStatsPlots\' , 'totalMean', '.tif'])
end

%%
% Now we cut the domain into sections of the size of an actin polymer in the simulation. Can find this info in the script used to run the simulation, hopefully in the same folder as this data. If it's not, go yell at whoever ran the simulation.

actin_length = 0.5;    % From script, variable "actin_length", length of actin monomers in units of um
nmonomer = 11;         % Number of monomers in each filament

polymerLength = floor(actin_length*nmonomer);  % Ensure that we use an integer. Unsure of whether to use floor or ceil. Can try both..
binedgesX = params.xRange(1):polymerLength:params.xRange(2); 
binedgesY = params.yRange(1):polymerLength:params.yRange(2);

binnedDivV = zeros(numel(binedgesX)-1, numel(binedgesY)-1, numel(time));

for ii=1:numel(binedgesX)-1
    xEdge = binedgesX(ii);
    indX = [grid.x >= xEdge & grid.x < xEdge+polymerLength];
    for jj=1:numel(binedgesY)-1
        local = divV; % Copy the original data to look at ony local info without destroying all original data
        yEdge = binedgesY(jj);
        indY = [grid.y >= yEdge & grid.y < yEdge+polymerLength];
        inThisBin = repmat(indX & indY, [1,1,numel(time)]);
        local(~inThisBin) = 0;
        binnedDivV(ii,jj,:) = mean(mean(local));
    end
end

save('interpedData', 'binnedDivV', '-append');

figure;
contractilityRatio = squeeze(min(min(binnedDivV)) ./ max(max(binnedDivV)));
stdRatio  =  std(contractilityRatio);
meanRatio = mean(contractilityRatio);
plot(time, contractilityRatio);
xlabel('time (s)');
ylabel('Min$\langle \vec{\nabla} \cdot \vec{v}_{local} \rangle$ / Max$\langle \vec{\nabla} \cdot \vec{v}_{local} \rangle$',...
    'Interpreter', 'latex')
legend(sprintf('\\mu \\pm \\sigma =  %.02e \\pm %.02e', meanRatio, stdRatio))
title({['Ratio of local divergences']; ['Bin Size = polymer length = ', num2str(polymerLength)]});
prettyFig;

if isunix
    saveas(gcf, [pwd, '/divStatsPlots/' , 'localRatio', '.tif'])
elseif ispc
    saveas(gcf, [pwd, '\divStatsPlots\' , 'localRatio', '.tif'])
end

%%
% Try smoothing the data
figure;
smoothRatio = smooth(contractilityRatio);
stdSmoothRatio  =  std(smoothRatio);
meanSmoothRatio = mean(smoothRatio);
plot(time, smoothRatio);
xlabel('time (s)');
ylabel('Min$\langle \vec{\nabla} \cdot \vec{v}_{local} \rangle$ / Max$\langle \vec{\nabla} \cdot \vec{v}_{local} \rangle$',...
    'Interpreter', 'latex')
legend(sprintf('\\mu \\pm \\sigma = %.02e \\pm %.02e', meanSmoothRatio, stdSmoothRatio))
title({['Ratio of local divergences (smoothed)']; ['Bin Size = polymer length = ', num2str(polymerLength)]});%, 'Interpreter','none')
prettyFig;

if isunix
    saveas(gcf, [pwd, '/divStatsPlots/' , 'localSmoothRatio', '.tif'])
elseif ispc
    saveas(gcf, [pwd, '\divStatsPlots\' , 'localSmoothRatio', '.tif'])
end

%%
% Now try looking only at the total divergence in an area in the center of the simulation domain

if isunix
    mkdir('divStatsPlots/reduceDomainSize');
elseif ispc
    mkdir('divStatsPlots\reduceDomainSize');
end

% Assuming a square matrix of divergence values, successively cut out first and last row/column, sum, plot, and repeat while domain gets smaller until only a (20 um)^2 region in the center

for ii=0:floor(size(divV,1)/2)-10
    innerDivV = divV(1+ii:end-ii, 1+ii:end-ii,:);
    innerMeanDivV = squeeze(mean(mean(innerDivV)));
    plot(time, innerMeanDivV)
    innerMeanTimeSeries = mean(innerMeanDivV);
    innerStdTimeSeries  =  std(innerMeanDivV);
    xlabel('time (s)')
    ylabel('$\langle \vec{\nabla} \cdot \vec{v} \rangle (t)$', 'Interpreter', 'latex')
    legend(sprintf('\\mu \\pm \\sigma = %.02e \\pm %.02e', innerMeanTimeSeries, innerStdTimeSeries), 'Location', 'Northwest')
    title({['Average over inner domain']; [sprintf('Center region (%d um)^2', params.L-2*ii*interpParams.dr)]});
    prettyFig;
    
    if isunix
        saveas(gcf, [pwd, '/divStatsPlots/reduceDomainSize/' , num2str(ii), '.tif'])
    elseif ispc
        saveas(gcf, [pwd, '\divStatsPlots\reduceDomainSize\' , num2str(ii), '.tif'])
    end
    clf;
end