clear; close all;
tic;
p = genpath('~/matlab/afinesAnalysis/');
addpath(p);


% First we just get all the data we can
cd('contractility')
aRange = 0:0.1:1;
pRange = 0:0.1:1;

[aGrid,pGrid] = meshgrid(aRange, pRange);

epsMat = zeros(size(aGrid));
epsdotMat = zeros(size(aGrid)); %initialize a couple of empty matrices to store phase space data
fileMat = cell(size(aGrid)); % also store files in same order so we know what's what...
maxOccupied = zeros(size(aGrid));
strainVars = {'maxStrain', 'maxStrainRate'};

for ii = 1:121
    aVal = aGrid(ii);
    pVal = pGrid(ii);
    fname = sprintf('a-%0.1f-p-%0.1f', aVal, pVal);
    cd(fname);
    pwd
    load('interpedData.mat')
    for jj=1:size(grid.n,3)
	occGrid(jj) = numel(find(grid.n(:,:,jj)));
    end
    epsMat(ii) = maxStrain;
    epsdotMat(ii) = maxStrainRate;
    fileMat{ii} = fname;
    cd ..;
end
cd ..
savedVars={'epsMat','epsdotMat', 'fileMat'};

if ispc
    save([pwd,'\phaseSpaceData_mWeighted.mat'], savedVars{:}); 
elseif isunix
    save([pwd,'/phaseSpaceData_mWeighted.mat'], savedVars{:});
end

% Plot colormap of strain
pcolor(aGrid,pGrid,epsMat);
xlabel('$Motor \; concentration (\mu m ^{-2})$', 'Interpreter','latex')
ylabel('$xLinker \; concentration (\mu m ^{-2})$', 'Interpreter', 'latex')
shading interp;
c = colorbar;
title(c, '$\epsilon_{max}$', 'Interpreter', 'latex')
saveas(gcf, 'strainPhaseSpace', 'fig');
saveas(gcf, 'strainPhaseSpace', 'tif');

% Plot colormap of strain rate
pcolor(aGrid,pGrid,epsdotMat);
xlabel('$Motor \; concentration (\mu m ^{-2})$', 'Interpreter','latex')
ylabel('$xLinker \; concentration (\mu m ^{-2})$', 'Interpreter', 'latex')
shading interp;
c = colorbar;
title(c, '$\dot{\epsilon}_{max} \; (s^{-1})$', 'Interpreter', 'latex')
saveas(gcf, 'strainRatePhaseSpace', 'fig');
saveas(gcf, 'strainRatePhaseSpace', 'tif');
