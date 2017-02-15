%%
% This script gets the statistics out of velocity divergence data output running divVelocity.m on a simdata.mat data file that was generated using read_data.m on the actins.txt output of an afines simulation
%
% Created by Daniel Seara at 2017/02/15 11:05

clear; close all;

fname='interpedData.mat';
myVars = {'grid', 'binParams', 'divV'};
load(fname, myVars{:});
fname='simdata.mat';
myVars={'adata','params'};
load(fname, myVars{:});

if ~exist('divStatsPlots','dir')
    mkdir('divStatsPlots')
end

if isunix
    file = strsplit(pwd,'/');
    file = file{end};
elseif ispc
    file = strsplit(pwd,'\');
    file = file{end};
end

% First we calculate the mean of the total velocity divergence over the whole domain and plot its time series

meanDivV = squeeze(mean(mean(divV)));
time = (params.timestep(1)+binParams.DeltaT*0.5):binParams.DeltaT:(params.timestep(end)-binParams.DeltaT*0.5); % Get the half-time steps where velocity is actually calculated
plot(time, meanDivV)
xlabel('time (s)')
ylabel('$\langle \vec{\nabla} \cdot \vec{v} \rangle (t)$', 'Interpreter', 'latex')
legend('Average over entire domain', 'Location', 'Northeast')
title(file, 'Interpreter','none')
if isunix
    saveas(gcf, [pwd, '/divStatsPlots/' , 'totalMean', '.tif'])
elseif ispc
    saveas(gcf, [pwd, '\divStatsPlots\' , 'totalMean', '.tif'])
end

%%
% Now we cut the domain into sections of the size of an actin polymer in the simulation. Can find this info in the script used to run the simulation, hopefully in the same folder as this data. If it's not, go yell at whoever ran the simulation.

monomerLength = 0.5;    % From script, length of actin monomers in units of um
nMonomer = 11;          % Number of monomers in each filament

polymerLength = 5;  % Ensure that we use an integer. Unsure of whether to use floor or ceil. Can try both..
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

contractilityRatio = squeeze(min(min(binnedDivV)) ./ max(max(binnedDivV)));
plot(time, contractilityRatio);
hold on;
plot(time, repmat(mean(contractilityRatio),[1,numel(time)]))