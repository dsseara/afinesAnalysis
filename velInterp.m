% Script for creating an interpolated velocity field from
% position of actin beads from AFiNeS simulation
% 
% Daniel Seara at 2017/01/09 21:54

clear;
close all;
tic

% Returns simdata.mat if doesn't exist already
%run([pwd, 'read_data2.m']);
fname = 'simdata.mat';
vars = {'adata', 'params'};
load(fname, vars{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set all parameters for interpolation
[numBeads, dim, numFrames] = size(adata);
binParams.numTimeSteps = 10; % How many frames between which to calculate velocity
binParams.DeltaT = binParams.numTimeSteps*params.dt;

% Bin velocities before interpolating to reduce noise
% bins of binSize and have a threshold number nThresh in each bin
binParams.binSize = 1; % in um
binParams.nThresh = 5; % Found by trial and error

interpParams.dr = binParams.binSize; % Size of grid to interpolate to
interpParams.interpRadius = binParams.binSize; % radius of interpolation region, anything beyond it is not considered
interpParams.d0 = 5 * binParams.binSize; % Freedman et al 2016, S2, width of Gaussian to weight (wider than interpRadius? Seems odd..)
interpParams.polygon = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xCoords = adata(:,1,1:binParams.numTimeSteps:numFrames);
yCoords = adata(:,2,1:binParams.numTimeSteps:numFrames);
xBases = squeeze(xCoords(:,:,1:end-1));
yBases = squeeze(yCoords(:,:,1:end-1));

d.x = squeeze(mod(diff(xCoords,1,3) + params.L/2, params.L) - params.L/2);
d.y = squeeze(mod(diff(yCoords,1,3) + params.L/2, params.L) - params.L/2);
%v.x = d.x./binParams.DeltaT;
%v.y = d.y./binParams.DeltaT;

% Make grid size mx2 of form [xg yg]m, grid size dr
[grid.x, grid.y] = meshgrid(params.xRange(1):interpParams.dr:params.xRange(2), params.yRange(1):interpParams.dr:params.yRange(2));
gridMat = [grid.x(:) grid.y(:)];

grid.dx = NaN([size(grid.x) size(d.x,2)]);
grid.dy = NaN([size(grid.y) size(d.y,2)]);

for i=1:size(d.x,2)
    disp(i)
    dxFrame = d.x(:,i);
    dyFrame = d.y(:,i);
    xbaseFrame = xBases(:,i);
    ybaseFrame = yBases(:,i);
    
    % Bin vectors within bins of size binSize, get average position and value if have at least nThresh beads in that bin
    dBinned  = binVectors(xbaseFrame,ybaseFrame,dxFrame,dyFrame,params.xRange,params.yRange,binParams.binSize,binParams.nThresh);

    % Interpolate frame
    dInterp = vectorFieldSparseInterpPatrick(dBinned, gridMat, interpParams.interpRadius, interpParams.d0, interpParams.polygon);
 
    interpedDX = reshape(dInterp(:,3), length(grid.x(1,:)), length(grid.x(:,1)));
    interpedDY = reshape(dInterp(:,4), length(grid.y(1,:)), length(grid.y(:,1)));
    interpedDX(isnan(interpedDX)) = 0;
    interpedDY(isnan(interpedDY)) = 0;

    grid.dx(:,:,i) = interpedDX;
    grid.dy(:,:,i) = interpedDY;
   
end % end loop over all the frames

clearvars -except interpParams binParams grid d
save([pwd,'/interpedData.mat'], '-append' )
toc