%%
% Script for calculating and the strain rate time series for Afines data
% This uses an external file "strain.m" which finds the strain rate tensor, eigenvectors
% and eigenvalues for the given vector field. 

clear; close all; tic;

fname = 'interpedData.mat';
myVars = {'grid','interpParams'};
load(fname, myVars{:});
fname = 'simdata.mat';
myVars = {'params'};
load(fname, myVars{:})

if isunix
    mkdir('Overlay/heatmap');
elseif ispc
    mkdir('Overlay\heatmap');
end

[width, height, nFrames] = size(grid.vx);
divV = zeros(size(grid.vx));
strain.tensor = zeros(width, height,2,2,nFrames);
strain.evecs = zeros(width, height,2,2,nFrames);
strain.evals = zeros(width, height,2,2,nFrames);


figure
for frame = 1:nFrames
    [tensor, evecs, evals] = symmetricGradient(grid.vx(:,:,frame), grid.vy(:,:,frame), interpParams.dr);
    strain.tensor(:,:,:,:,frame) = tensor;
    strain.evecs(:,:,:,:,frame) = evecs;
    strain.evals(:,:,:,:,frame) = evals;
    divDR(:,:,frame) = sum(sum(evals,3),4);

    % pcolor(grid.x, grid.y, divDR(:,:,frame))
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

    % if isunix
    %     saveas(gcf, [pwd, '/Overlay/heatmap/' ,num2str(frame), '.tif'])
    % elseif ispc
    %     saveas(gcf, [pwd, '\Overlay\heatmap\' ,num2str(frame), '.tif'])
    % end

    % clf
end

savedVars = {'divDR', 'strain'};
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