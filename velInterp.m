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
vars = {'adata', 'timestep'};
load(fname, vars{:});

% Change this if you want to find velocities for steps longer than 
% the time between frames
numTimeSteps = 1;

% Enforce periodic boundary conditions
xyDisplacement = mod(diff(adata(:,1:2, 1:numTimeSteps:end),1,3) + L/2, L) - L/2;

[numBeads, dim, numFrames] = size(xyDisplacement);
dt = numTimeSteps*uniquetol(diff(timestep));

xCoords = adata(:,1,1:numTimeSteps:end);
yCoords = adata(:,2,1:numTimeSteps:end); clear adata;
xBases = squeeze(xCoords(:,:,1:end-1)); clear xCoords;
yBases = squeeze(yCoords(:,:,1:end-1)); clear yCoords;
xRange = [floor(min(min(xBases))) ceil(max(max(xBases)))];
yRange = [floor(min(min(yBases))) ceil(max(max(yBases)))];

dx = squeeze(xyDisplacement(:,1,:));
dy = squeeze(xyDisplacement(:,2,:)); clear xyDisplacement;
vx = dx./dt;
vy = dy./dt;

% Bin velocities before interpolating to reduce noise
% bins of binSize and have a threshold number nThresh in each bin
binSize=1; % in um
nThresh = 5; % Found by trial and error

% Make grid size mx2 of form [xg yg]m, grid size dr
dr = binSize;
[xGrid, yGrid] = meshgrid(xRange(1):dr:xRange(2), yRange(1):dr:yRange(2));
gridMat = [xGrid(:) yGrid(:)];


interpRadius = 2*binSize; % radius of interpolation region, anything beyond it is not considered
d0 = 5*binSize; % Freedman et al 2016, S2, width of Gaussian to weight (wider than interpRadius? Seems odd..)
polygon = [];

adataVelInterped = NaN([size(xGrid) dim numFrames]);

for i=1:numFrames
    disp(i)
    vxFrame = vx(:,i);
    vyFrame = vy(:,i);
    xbaseFrame = xBases(:,i);
    ybaseFrame = yBases(:,i);
    
    % Bin vectors within bins of size binSize, get average position and value if have at least nThresh beads in that bin
    vBinned  = binVectors(xbaseFrame,ybaseFrame,vxFrame,vyFrame,xRange,yRange,binSize,nThresh);

    % Interpolate frame
    vInterp = vectorFieldSparseInterpPatrick(vBinned,  gridMat, interpRadius, d0, polygon);
 
    interpedVX = reshape(vInterp(:,3), length(xGrid(1,:)), length(xGrid(:,1)));
    interpedVY = reshape(vInterp(:,4), length(yGrid(1,:)), length(yGrid(:,1)));
    interpedVX(isnan(interpedVX)) = 0;
    interpedVY(isnan(interpedVY)) = 0;

    adataVelInterped(:,:,1,i) = interpedVX;
    adataVelInterped(:,:,2,i) = interpedVY;
   
end % end loop over all the frames
save([pwd,'/interpedData.mat'] )
toc