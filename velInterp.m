% Script for creating an interpolated velocity field from
% position of actin beads from AFiNeS simulation

% Daniel Seara at 2017/01/09 21:54

clear;
close all;

% Returns simdata.mat if doesn't exist already
%run([pwd, 'read_data2.m']);
fname = 'simdata.mat';
vars = {'adata', 'timestep'};
load(fname, vars{:});


% Change this if you want to find velocities for steps longer than 
% the time between frames
numTimeSteps = 1;

xyDisplacement = diff(adata(:,1:2, 1:numTimeSteps:end),1,3);
[numBeads, dim, numFrames] = size(xyDisplacement);
dt = uniquetol(diff(timestep));
dt = repmat(dt, size(xyDisplacement));

adataVel = xyDisplacement ./ dt;

%% Bin all the vectors according to their starting positions in boxes of size L, [L] = μm
%  replace with average starting position and velocity to reduce noise

L=1;

xCoords = adata(:,1,1:numTimeSteps:end);
yCoords = adata(:,2,1:numTimeSteps:end);
vX = adataVel(:,1,:);
vY = adataVel(:,2,:);

%clear adata, adataVel

minX = floor(min(min(xCoords)));
minY = floor(min(min(yCoords)));
maxX =  ceil(max(max(xCoords)));
maxY =  ceil(max(max(yCoords)));

xRange = maxX - minX;
yRange = maxY - minY;

%bincenters = 



% %% Ignore 10% of linear range around the edges to avoid effects from periodic boundary conditions
% indX = xCoords > (minX + 0.05*xRange) & xCoords < (maxX - 0.05*xRange);
% indY = yCoords > (minY + 0.05*yRange) & yCoords < (maxY - 0.05*yRange);
% %%

% Create vector matrix to be input into vectorFieldSparseInterpPatrick.m
% nx4 matrix of form [y0 x0 y x]n, y0 x0 are base of vectors, y,x are tips
% base is given by xCoords and yCoords, tips by xCoords+vx, similarly for y

adataVelFormatted = [yCoords(:,:,1:end-1) xCoords(:,:,1:end-1) yCoords(:,:,1:end-1)+vY xCoords(:,:,1:end-1)+vX];

% Make grid size mx2 of form [yg xg]m, grid size dr

dr = 1;
[xGrid, yGrid] = meshgrid(minX:dr:maxX, minY:dr:maxY);
gridMat = [yGrid(:) xGrid(:)];

threshold = 1;
d0 = 5*dr; % Freedman et al 2016, S2
polygon = [];

% Interpolate test frame
frame = floor(numFrames*0.5);
testVec = adataVelFormatted(:,:,frame);
testInterp = vectorFieldSparseInterpPatrick(testVec, gridMat, threshold, d0, polygon);

interpedVX = reshape(testInterp(:,4), length(minX:dr:maxX), length(minX:dr:maxX));
interpedVY = reshape(testInterp(:,3), length(minY:dr:maxY), length(minY:dr:maxY));

% Now we compare the results of the interpolated field with the original one
figure()
quiver(xGrid, yGrid, interpedVX, interpedVY,5)
title('Interpolated')
figure()
quiver(squeeze(xCoords(:,:,frame)), squeeze(yCoords(:,:,frame)), squeeze(vX(:,:,frame)), squeeze(vY(:,:,frame)),5)
title('Original')
% legend('Interpolated', 'Original')
% title('Velocity Field interpolation test')
% hold off