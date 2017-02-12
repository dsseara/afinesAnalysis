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

dx = squeeze(mod(diff(xCoords,1,3) + params.L/2, params.L) - params.L/2);
dy = squeeze(mod(diff(yCoords,1,3) + params.L/2, params.L) - params.L/2);
v.x = dx./binParams.DeltaT;
v.y = dy./binParams.DeltaT;

% Make grid size mx2 of form [xg yg]m, grid size dr
[grid.x, grid.y] = meshgrid(params.xRange(1):interpParams.dr:params.xRange(2), params.yRange(1):interpParams.dr:params.yRange(2));
gridMat = [grid.x(:) grid.y(:)];

grid.vx = NaN([size(grid.x) size(v.x,2)]);
grid.vy = NaN([size(grid.y) size(v.y,2)]);

for i=1:size(v.x,2)
    disp(i)
    vxFrame = v.x(:,i);
    vyFrame = v.y(:,i);
    xbaseFrame = xBases(:,i);
    ybaseFrame = yBases(:,i);
    
    % Bin vectors within bins of size binSize, get average position and value if have at least nThresh beads in that bin
    vBinned  = binVectors(xbaseFrame,ybaseFrame,vxFrame,vyFrame,params.xRange,params.yRange,binParams.binSize,binParams.nThresh);

    % Interpolate frame
    vInterp = vectorFieldSparseInterpPatrick(vBinned, gridMat, interpParams.interpRadius, interpParams.d0, interpParams.polygon);
 
    interpedVX = reshape(vInterp(:,3), length(grid.x(1,:)), length(grid.x(:,1)));
    interpedVY = reshape(vInterp(:,4), length(grid.y(1,:)), length(grid.y(:,1)));
    interpedVX(isnan(interpedVX)) = 0;
    interpedVY(isnan(interpedVY)) = 0;

    grid.vx(:,:,i) = interpedVX;
    grid.vy(:,:,i) = interpedVY;
   
end % end loop over all the frames

clearvars -except interpParams binParams grid v
save([pwd,'/interpedData3.mat'] )
toc