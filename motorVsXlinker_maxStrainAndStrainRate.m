% This script analyzes afines output data to investigate the max strain and strain rate
% in a contractile actomyosin gel as a function of varying motor and cross-linker concentrations.
% Assumes that there are a list of files with names of the form a_x.x_p_y.y, where x.x and y.y 
% are the concentrations used in the enclosing folder, from 0.0-1.0 in 0.1 increments.
%
% Created by Daniel Seara at 2017/01/30 15:47

%%
% First we just get all the data we can
tic;
aRange = 0:0.1:1;
pRange = 0:0.1:1;

[aGrid,pGrid] = meshgrid(aRange, pRange);

epsMat = zeros(size(aGrid));
epsdotMat = zeros(size(aGrid)); %initialize a couple of empty matrices to store phase space data
fileMat = cell(size(aGrid)); % also store files in same order so we know what's what...

ndt = 1; % This controls what the number of time steps are between calculated positions

for ii = 1:numel(aGrid)
    aVal = aGrid(ii);
    pVal = pGrid(ii);
    fname = sprintf('a_%0.1f_p_%0.1f', aVal, pVal);
    [epsmax, epsdotmax] = strainFromDisplacement(fname, ndt);
    epsMat(ii) = epsmax;
    epsdotMat(ii) = epsdotmax;
    fileMat{ii} = fname;
end

clearvars -except epsMat epsdotMat fileMat

if ispc
    save([pwd,'\phaseSpaceData.mat'], 'epsMat', 'epsdotMat', 'fileMat');
elseif isunix
    save([pwd,'/phaseSpaceData.mat'], 'epsMat', 'epsdotMat', 'fileMat');
end

clear;
toc

% Copied some code prior to clear statement for plotting (IAL: 1/31/17)
aRange = 0:0.1:1;
pRange = 0:0.1:1;
[aGrid,pGrid] = meshgrid(aRange, pRange);

load phaseSpaceData.mat

% Plot colormap of strain
pcolor(aGrid,pGrid,epsMat);
xlabel('$Active \; motor \; concentration (\mu m ^{-2})$', 'Interpreter','latex')
ylabel('$Passive \; motor \; concentration (\mu m ^{-2})$', 'Interpreter', 'latex')
shading interp;
c = colorbar;
ylabel(c, '$\epsilon_{max}$', 'Interpreter', 'latex')
saveas(gcf, 'strainPhaseSpace', 'fig');
saveas(gcf, 'strainPhaseSpace', 'epsc');

% Plot colormap of strain rate
pcolor(aGrid,pGrid,epsdotMat);
xlabel('$Active \; motor \; concentration (\mu m ^{-2})$', 'Interpreter','latex')
ylabel('$Passive \; motor \; concentration (\mu m ^{-2})$', 'Interpreter', 'latex')
shading interp;
c = colorbar;
ylabel(c, '$\dot{\epsilon}_{max} \; (s^{-1})$', 'Interpreter', 'latex')
saveas(gcf, 'strainRatePhaseSpace', 'fig');
saveas(gcf, 'strainRatePhaseSpace', 'epsc');


