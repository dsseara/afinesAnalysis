function [maxStrain, strainRate2575] = divStats()
% This script gets the statistics out of velocity divergence data output
% running divVelocity.m on a simdata.mat data file that was generated using
% read_data.m on the actins.txt output of an afines simulation
%
%
% Created by Daniel Seara at 2017/02/15 11:05

clear; close all;

fname  = 'interpedData.mat';
myVars = {'grid', 'binParams', 'divV', 'divDR', 'interpParams'};
load(fname, myVars{:});

fname  = 'simdata.mat';
myVars = {'params'};
load(fname, myVars{:});

time = params.timestep(1:end-binParams.numTimeSteps); 

% Get size of polymers as a length scale...
actin_length  = 0.5;    % From script, variable "actin_length", length of actin monomers in units of um
nmonomer      = 11;     % Number of monomers in each filament
polymerLength = floor(actin_length*nmonomer);  % Ensure that we use an integer. Unsure of whether to use floor or ceil. Can try both..

% Look at everything between [-25, 25] in both directions
innerXInd = grid.x(1,:)>=(-25) & grid.x(1,:)<=(25);
innerYInd = grid.y(:,1)>=(-25) & grid.y(:,1)<=(25);

innerStrain     = cumStrain(innerXInd, innerYInd, :);
innerMeanStrain = squeeze(mean(mean(innerStrain)));

% Get some statistics about the strain, the max and the slope to go from 25%-75% of the max value
maxStrain = max(innerMeanStrain);
[~, ind25] = min(abs(innerMeanStrain - 0.25*maxStrainAbs));
[~, ind75] = min(abs(innerMeanStrain - 0.75*maxStrainAbs));
strainRate2575 = (innerMeanStrain(ind75) - innerMeanStrain(ind25))./(time(ind75) - time(ind25));

% Plot it up
plot(time, innerMeanStrain)
xlabel('time (s)')
ylabel('$\langle |\epsilon|\rangle (t)$', 'Interpreter', 'latex')
legend(sprintf('Max = %.02f, .25-.75 slope = %.02f', maxStrain, strainRate2575), 'Location', 'SouthEast')
title({['Average strain']; ['Center region (50 um)^2']});
prettyFig;

if isunix
    saveas(gcf, [pwd, '/strain.tif'])
    saveas(gcf, [pwd, '/strain.fig'])
elseif ispc
    saveas(gcf, [pwd, '\strain.tif'])
    saveas(gcf, [pwd, '\strain.fig'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% This will always give zero, don't do it! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Because of mass conservation, div V = 0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % First we calculate the mean of the total velocity divergence over the whole domain and plot its time series
% 
% meanDivV  = squeeze(mean(mean(divV)));
% plot(time, meanDivV)
% meanTimeSeries = mean(meanDivV);
% stdTimeSeries  =  std(meanDivV);
% xlabel('time (s)')
% ylabel('$\langle \vec{\nabla} \cdot \vec{v} \rangle (t)$', 'Interpreter', 'latex')
% legend(sprintf('\\mu \\pm \\sigma = %.02e \\pm %.02e', meanTimeSeries, stdTimeSeries))
% title('Average over whole domain')%, 'Interpreter','none')
% prettyFig;
% 
% if isunix
%     saveas(gcf, [pwd, '/divStatsPlots/' , 'strainRate_totalMean', '.tif'])
% elseif ispc
%     saveas(gcf, [pwd, '\divStatsPlots\' , 'strainRate_totalMean', '.tif'])
% end
% 
% 
% % Now we do the same thing for the strain
% meanStrain  = squeeze(mean(mean(cumStrain)));
% plot(time, meanStrain)
% meanTimeSeries = mean(meanStrain);
% stdTimeSeries  =  std(meanStrain);
% xlabel('time (s)')
% ylabel('$\langle | \epsilon | \rangle (t)$', 'Interpreter', 'latex')
% legend(sprintf('\\mu \\pm \\sigma = %.02e \\pm %.02e', meanTimeSeries, stdTimeSeries))
% title('Average over whole domain')%, 'Interpreter','none')
% prettyFig;
% 
% if isunix
%     saveas(gcf, [pwd, '/divStatsPlots/' , 'strain_totalMean', '.tif'])
% elseif ispc
%     saveas(gcf, [pwd, '\divStatsPlots\' , 'strain_totalMean', '.tif'])
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Now we cut the domain into sections of the size of an actin polymer in the simulation. 
% % Can find this info in the script used to run the simulation, hopefully in the same folder as this data. 
% % If it's not, go yell at whoever ran the simulation.
% 
% binedgesX = params.xRange(1):polymerLength:params.xRange(2); 
% binedgesY = params.yRange(1):polymerLength:params.yRange(2);
% 
% binnedDivV = zeros(numel(binedgesX)-1, numel(binedgesY)-1, numel(time));
% 
% for ii=1:numel(binedgesX)-1
%     xEdge = binedgesX(ii);
%     indX = [grid.x >= xEdge & grid.x < xEdge+polymerLength];
%     for jj=1:numel(binedgesY)-1
%         local = divV; % Copy the original data to look at ony local info without destroying all original data
%         yEdge = binedgesY(jj);
%         indY = [grid.y >= yEdge & grid.y < yEdge+polymerLength];
%         inThisBin = repmat(indX & indY, [1,1,numel(time)]);
%         local(~inThisBin) = 0;
%         binnedDivV(ii,jj,:) = mean(mean(local));
%     end
% end
% 
% save('interpedData', 'binnedDivV', '-append');
% 
% figure;
% contractilityRatio = squeeze(min(min(binnedDivV)) ./ max(max(binnedDivV)));
% stdRatio  =  std(contractilityRatio);
% meanRatio = mean(contractilityRatio);
% plot(time, contractilityRatio);
% xlabel('time (s)');
% ylabel('Min$\langle \vec{\nabla} \cdot \vec{v}_{local} \rangle$ / Max$\langle \vec{\nabla} \cdot \vec{v}_{local} \rangle$',...
%     'Interpreter', 'latex')
% legend(sprintf('\\mu \\pm \\sigma =  %.02e \\pm %.02e', meanRatio, stdRatio))
% title({['Ratio of local divergences']; ['Bin Size = polymer length = ', num2str(polymerLength)]});
% prettyFig;
% 
% if isunix
%     saveas(gcf, [pwd, '/divStatsPlots/' , 'localRatio', '.tif'])
% elseif ispc
%     saveas(gcf, [pwd, '\divStatsPlots\' , 'localRatio', '.tif'])
% end
% 
% %%
% % Try smoothing the data
% figure;
% smoothRatio = smooth(contractilityRatio);
% stdSmoothRatio  =  std(smoothRatio);
% meanSmoothRatio = mean(smoothRatio);
% plot(time, smoothRatio);
% xlabel('time (s)');
% ylabel('Min$\langle \vec{\nabla} \cdot \vec{v}_{local} \rangle$ / Max$\langle \vec{\nabla} \cdot \vec{v}_{local} \rangle$',...
%     'Interpreter', 'latex')
% legend(sprintf('\\mu \\pm \\sigma = %.02e \\pm %.02e', meanSmoothRatio, stdSmoothRatio))
% title({['Ratio of local divergences (smoothed)']; ['Bin Size = polymer length = ', num2str(polymerLength)]});%, 'Interpreter','none')
% prettyFig;
% 
% if isunix
%     saveas(gcf, [pwd, '/divStatsPlots/' , 'localSmoothRatio', '.tif'])
% elseif ispc
%     saveas(gcf, [pwd, '\divStatsPlots\' , 'localSmoothRatio', '.tif'])
% end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Now try looking only at the total divergence in an area in the center of the simulation domain

% if isunix
%     mkdir('divStatsPlots/strain');
%     mkdir('divStatsPlots/strainrate');
% elseif ispc
%     mkdir('divStatsPlots\strain');
%     mkdir('divStatsPlots\strainrate');
% end

% Assuming a square matrix of divergence values, successively cut out first and last row/column, sum, plot, and repeat while domain gets smaller until only a (20 um)^2 region in the center

% Look at the center (50um)^2 region of the (100um)^2 total domain, i.e. cut out 25 um of 

% Look at progressively smaller regions
% Whole domain is [-L/2, L/2], so start there and cut smaller and smaller to [-L/10, L/10]
% fracs = 2:0.5:10;

% for ii=fracs

%     L_temp = (params.L/ii)*2;
%     innerXInd = grid.x(1,:)>=(-params.L/ii) & grid.x(1,:)<=(params.L/ii);
%     innerYInd = grid.y(:,1)>=(-params.L/ii) & grid.y(:,1)<=(params.L/ii);

%     innerDivV = divV(innerXInd, innerYInd, :);
%     innerMeanDivV = squeeze(mean(mean(innerDivV)));
%     plot(time, innerMeanDivV)
%     innerMeanTimeSeries = mean(innerMeanDivV);
%     innerStdTimeSeries  =  std(innerMeanDivV);
%     xlabel('time (s)')
%     ylabel('$\langle \vec{\nabla} \cdot \vec{v} \rangle (t)$', 'Interpreter', 'latex')
%     legend(sprintf('\\mu \\pm \\sigma = %.02e \\pm %.02e', innerMeanTimeSeries, innerStdTimeSeries), 'Location', 'Northwest')
%     title({['Average strain rate']; [sprintf('Center region (%0.02f um)^2', L_temp)]});
%     prettyFig;

%     %save('interpedData','innerDivV', '-append')

%     if isunix
%         saveas(gcf, [pwd, '/divStatsPlots/strainrate/' , num2str(L_temp), '_strainRate.tif'])
%     elseif ispc
%         saveas(gcf, [pwd, '\divStatsPlots\strainrate\' , num2str(L_temp), '_strainRate.tif'])
%     end
%     clf;

% % Same thing for strain

%     innerStrain = cumStrain(innerXInd, innerYInd, :);
%     innerMeanStrain = squeeze(mean(mean(innerStrain)));
%     plot(time, innerMeanStrain)
%     % Get some statistics about the strain, the max and the slope to go from 25%-75% of the max value
%     maxStrain = max(innerMeanStrain);
    
%     %innerMeanTimeSeries = mean(innerMeanStrain);
%     %innerStdTimeSeries  =  std(innerMeanStrain);
%     xlabel('time (s)')
%     ylabel('$\langle |\epsilon|\rangle (t)$', 'Interpreter', 'latex')
%     %legend(sprintf('\\mu \\pm \\sigma = %.02e \\pm %.02e', innerMeanTimeSeries, innerStdTimeSeries), 'Location', 'Northwest')
%     title({['Average strain']; [sprintf('Center region (%0.02f um)^2', L_temp)]});
%     prettyFig;
    

%     %save('interpedData','innerStrain', '-append')

%     if isunix
%         saveas(gcf, [pwd, '/divStatsPlots/strain/' , num2str(L_temp), '_strain.tif'])
%     elseif ispc
%         saveas(gcf, [pwd, '\divStatsPlots\strain\' , num2str(L_temp), '_strain.tif'])
%     end
%     clf;
% end