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

% Change this if you want to find velocities for steps longer than 
% the time between frames
interpParams.dt = 10;

xCoords = adata(:,1,1:numTimeSteps:end);
yCoords = adata(:,2,1:numTimeSteps:end);
xBases = squeeze(xCoords(:,:,1:end-1)); clear xCoords;
yBases = squeeze(yCoords(:,:,1:end-1)); clear yCoords;

xyDisplacement = mod(diff(adata(:,1:2, 1:numTimeSteps:end),1,3) + params.L/2, params.L) - params.L/2; clear adata;

[numBeads, dim, numFrames] = size(xyDisplacement);

dx = squeeze(xyDisplacement(:,1,:));
dy = squeeze(xyDisplacement(:,2,:)); clear xyDisplacement;
v.x = dx./interpParams.dt; clear dx;
v.y = dy./interpParams.dt; clear dy;

% Bin velocities before interpolating to reduce noise
% bins of binSize and have a threshold number nThresh in each bin
binSize=1; % in um
nThresh = 5; % Found by trial and error

% Make grid size mx2 of form [xg yg]m, grid size dr
dr = binSize;
[xGrid, yGrid] = meshgrid(params.xRange(1):dr:params.xRange(2), params.yRange(1):dr:params.yRange(2));
gridMat = [xGrid(:) yGrid(:)];


interpRadius = binSize; % radius of interpolation region, anything beyond it is not considered
d0 = 5*binSize; % Freedman et al 2016, S2, width of Gaussian to weight (wider than interpRadius? Seems odd..)
polygon = [];

adataVelInterped = NaN([size(xGrid) dim numFrames]);

for i=1:numFrames
    disp(i)
    vxFrame = v.x(:,i);
    vyFrame = v.y(:,i);
    xbaseFrame = xBases(:,i);
    ybaseFrame = yBases(:,i);
    
    % Bin vectors within bins of size binSize, get average position and value if have at least nThresh beads in that bin
    vBinned  = binVectors(xbaseFrame,ybaseFrame,vxFrame,vyFrame,xRange,yRange,binSize,nThresh);

    % Interpolate frame
    vInterp = vectorFieldSparseInterpPatrick(vBinned, gridMat, interpRadius, d0, polygon);
 
    interpedVX = reshape(vInterp(:,3), length(xGrid(1,:)), length(xGrid(:,1)));
    interpedVY = reshape(vInterp(:,4), length(yGrid(1,:)), length(yGrid(:,1)));
    interpedVX(isnan(interpedVX)) = 0;
    interpedVY(isnan(interpedVY)) = 0;

    adataVelInterped(:,:,1,i) = interpedVX;
    adataVelInterped(:,:,2,i) = interpedVY;
   
end % end loop over all the frames
save([pwd,'/interpedData2.mat'] )
toc