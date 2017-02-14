%%
% Script for calculating and the strain rate time series for Afines data
% This uses an external file "strain.m" which finds the strain rate tensor, eigenvectors
% and eigenvalues for the given vector field. 

clear; close all; tic;

fname = 'interpedData.mat';
myVars = {'grid','interpParams'};
load(fname, myVars{:});
fname = 'simdata.mat';
myVars = {'adata', 'params'};
load(fname, myVars{:})


[~, ~, nFrames] = size(grid.vx);
divV = zeros(size(grid.vx));
strainRate.tensor = zeros(51,51,2,2,nFrames);
strainRate.evecs = zeros(51,51,2,2,nFrames);
strainRate.evals = zeros(51,51,2,2,nFrames);

% We are going to want to cut the system into an NxN grid and calculate the local divergence in each section. To do this we will first make sure that the given N actually cuts the grid up nicely
N = 5;
if mod(params.L/interpParams.dr,N) ~= 0
    error('Currently trying to make %dx%d grid, but number of grid points aren"t divisible by %d', N,N,N)
end
figure
for frame = 1:nFrames
    [tensor, evecs, evals] = strain(grid.vx(:,:,frame), grid.vy(:,:,frame), interpParams.dr);
    strainRate.tensor(:,:,:,:,frame) = tensor;
    strainRate.evecs(:,:,:,:,frame) = evecs;
    strainRate.evals(:,:,:,:,frame) = evals;
    divV(:,:,frame) = sum(sum(evals,3),4);

    % pcolor(grid.x, grid.y, divV(:,:,frame))
    % shading interp
    % colormap jet
    % caxis([-1.5, 1.5])
    % colorbar;
    % hold on;
    % quiver(grid.x, grid.y, grid.vx(:,:,frame), grid.vy(:,:,frame),'k')
    
    % xlim([params.xRange(1) params.xRange(2)])
    % ylim([params.yRange(1) params.yRange(2)])
    % pbaspect([1 1 1])
    % set(gca,'xtick',[])
    % set(gca,'xticklabel',[])
    % set(gca,'ytick',[])
    % set(gca,'yticklabel',[])
    % axis tight 
    % saveas(gcf, [pwd, '/Overlay/heatmap/' ,num2str(frame), '.tif'])
    % clf
end

savedVars = {'divV', 'strainRate'};
save('interpedData',savedVars{:}, '-append')


% errorbar(0:(numel(seriesStats.mean)-1), seriesStats.mean, seriesStats.std, 'o');
% hold on;
% plot([0,30],[0,0])
% xlim([0,30])
% legend('Mean Strain Rate, \mu \pm \sigma', 'Zero Line')
% seriesStats.std = [seriesStats.std std(strainRateSeries)];
% seriesStats.mean = [seriesStats.mean mean(strainRateSeries)];
% plot((1:10:2000)./4,strainRateSeries(1:10:end),'o')
% xlabel('time (s)')
% ylabel('$ \dot{\epsilon}(t) $', 'Interpreter', 'LaTex')
% txt = ['Mean $\pm$ std = ', num2str(seriesStats.mean(end)),'$\pm$' num2str(seriesStats.std(end))];
% text(250, 0.5*max(strainRateSeries), txt, 'Interpreter','LaTeX');
% txt2 = ['X-range = ', num2str(interiorXRange(1)), '-', num2str(interiorXRange(2))];
% txt3 = [', Y-range = ', num2str(interiorYRange(1)), '-', num2str(interiorYRange(2))];
% title([txt2, txt3])