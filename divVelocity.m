%%
% Script for calculating and the divergence of the velocity and displacement fields for Afines data
% 
% External files:   simdata.mat
%                   interpedData.mat
%                   symmetricGradient.m
%
% Created by Daniel Seara at 2017/02/16 12:19

clear; close all; tic;

fname = 'interpedData.mat';
myVars = {'grid','interpParams'};
load(fname, myVars{:});

[width, height, nFrames] = size(grid.vx);
divV  = zeros(size(grid.vx));
divDR = zeros(size(grid.dx));

strain.tensor     = zeros(width, height,2,2,nFrames);
strain.evecs      = zeros(width, height,2,2,nFrames);
strain.evals      = zeros(width, height,2,2,nFrames);
strainRate.tensor = zeros(width, height,2,2,nFrames);
strainRate.evecs  = zeros(width, height,2,2,nFrames);
strainRate.evals  = zeros(width, height,2,2,nFrames);

for frame = 1:nFrames
    [tensor, evecs, evals] = symmetricGradient(grid.dx(:,:,frame), grid.dy(:,:,frame), interpParams.dr);
    strain.tensor(:,:,:,:,frame) = tensor;
    strain.evecs(:,:,:,:,frame)  = evecs;
    strain.evals(:,:,:,:,frame)  = evals;
    divDR(:,:,frame) = sum(sum(evals,3),4);

    [tensor, evecs, evals] = symmetricGradient(grid.vx(:,:,frame), grid.vy(:,:,frame), interpParams.dr);
    strainRate.tensor(:,:,:,:,frame) = tensor;
    strainRate.evecs(:,:,:,:,frame)  = evecs;
    strainRate.evals(:,:,:,:,frame)  = evals;
    divV(:,:,frame) = sum(sum(evals,3),4);
end

savedVars = {'divDR', 'strain', 'divV', 'strainRate'};
save('interpedData', savedVars{:}, '-append')