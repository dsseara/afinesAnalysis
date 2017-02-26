%%
% Script that interpolates a vector field onto a grid using a Gaussian RBF Function for
% AFiNeS simulation outputs
%
% Created by Daniel Seara at 2017/02/26 13:36

clear;
close all;
tic

% Returns simdata.mat if doesn't exist already
run run_read_data.m;

fname = 'simdata.mat';
vars = {'adata', 'params'};
load(fname, vars{:});
[numBeads, dim, numFrames] = size(adata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bin velocities before interpolating to reduce noise
% bins of binSize and have a threshold number nThresh in each bin
binParams.numTimeSteps = 10; % How many frames between which to calculate velocity
binParams.DeltaT = binParams.numTimeSteps*params.dt;
binParams.binSize = 1; % in um
binParams.nThresh = 10; % Found by trial and error

% Set all parameters for interpolation
interpParams.dr = 2; % Size of grid to interpolate to in um
interpParams.interpRadius = interpParams.dr; % radius of interpolation region, anything beyond it is not considered
interpParams.d0 = 5 * binParams.binSize; % Freedman et al 2016, S2, width of Gaussian to weight (wider than interpRadius? Seems odd..)
interpParams.polygon = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xCoords = squeeze(adata(:,1,:));
yCoords = squeeze(adata(:,2,:));

dx = xCoords(:,1+binParams.numTimeSteps:end) - xCoords(:, 1:end-binParams.numTimeSteps);
dy = yCoords(:,1+binParams.numTimeSteps:end) - yCoords(:, 1:end-binParams.numTimeSteps);

% Enforce periodicity

d.x = mod(dx + params.L/2, params.L) - params.L/2;
d.y = mod(dy + params.L/2, params.L) - params.L/2;
v.x = d.x./binParams.DeltaT;
v.y = d.y./binParams.DeltaT;

% Make grid size mx2 of form [xg yg]m, grid size dr
[grid.x, grid.y] = meshgrid(params.xRange(1):interpParams.dr:params.xRange(2), params.yRange(1):interpParams.dr:params.yRange(2));
gridMat = [grid.x(:) grid.y(:)];

dx = NaN([size(grid.x) size(d.x,2)]);
dy = NaN([size(grid.y) size(d.y,2)]);
vx = NaN([size(grid.x) size(d.x,2)]);
vy = NaN([size(grid.y) size(d.y,2)]);
[height, width] = size(grid.x);

parfor frame=1:size(d.x,2)
    disp(['frame = ', num2str(frame)])
    dxFrame = d.x(:,frame);
    dyFrame = d.y(:,frame);
    vxFrame = v.x(:,frame);
    vyFrame = v.y(:,frame);
    xbaseFrame = xCoords(:,frame);
    ybaseFrame = yCoords(:,frame);
    
    % Bin vectors within bins of size binSize, get average position and value if have at least nThresh beads in that bin
    dBinned = binVectors(xbaseFrame, ybaseFrame, dxFrame, dyFrame, params.xRange, params.yRange, binParams.binSize, binParams.nThresh);
    vBinned = binVectors(xbaseFrame, ybaseFrame, vxFrame, vyFrame, params.xRange, params.yRange, binParams.binSize, binParams.nThresh);

    % Interpolate frame
    dInterp = vectorFieldSparseInterpPatrick(dBinned, gridMat, interpParams.interpRadius, interpParams.d0, interpParams.polygon);
    vInterp = vectorFieldSparseInterpPatrick(vBinned, gridMat, interpParams.interpRadius, interpParams.d0, interpParams.polygon);
 
    interpedDX = reshape(dInterp(:,3), height, width);
    interpedDY = reshape(dInterp(:,4), height, width);
    interpedDX(isnan(interpedDX)) = 0;
    interpedDY(isnan(interpedDY)) = 0;

    interpedVX = reshape(vInterp(:,3), height, width);
    interpedVY = reshape(vInterp(:,4), height, width);
    interpedVX(isnan(interpedVX)) = 0;
    interpedVY(isnan(interpedVY)) = 0;

    dx(:,:,frame) = interpedDX;
    dy(:,:,frame) = interpedDY;
    vx(:,:,frame) = interpedVX;
    vy(:,:,frame) = interpedVY;
end % end parfor loop over all the frames

grid.dx = dx;
grid.dy = dy;
grid.vx = vx;
grid.vy = vy;

clearvars -except grid d v binParams interpParams

save([pwd,'interpedData.mat'])
toc