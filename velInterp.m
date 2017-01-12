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

xyDisplacement = diff(adata(:,1:2, 1:numTimeSteps:end),1,3);
[numBeads, dim, numFrames] = size(xyDisplacement);
dt = uniquetol(diff(timestep));
dt = repmat(dt, size(xyDisplacement));

adataVel = xyDisplacement ./ dt;

%% Bin all the vectors according to their starting positions in boxes of size L, [L] = μm
%  replace with average starting position and velocity to reduce noise

xCoords = adata(:,1,1:numTimeSteps:end);
yCoords = adata(:,2,1:numTimeSteps:end);
xBases = squeeze(xCoords(:,:,1:end-1));
yBases = squeeze(yCoords(:,:,1:end-1));
vX = squeeze(adataVel(:,1,:));
vY = squeeze(adataVel(:,2,:));

%clear adata adataVel;

minX = floor(min(min(xCoords)));
minY = floor(min(min(yCoords)));
maxX =  ceil(max(max(xCoords)));
maxY =  ceil(max(max(yCoords)));

% Test of bin velocity function
% bins of binSize and have a threshold number nThresh in each bin
binSize=1;
nThresh = 0;

testFrame = floor(0.1*numFrames);
vxTest = vX(:,testFrame);
vyTest = vY(:,testFrame);
xbaseTest = xBases(:,testFrame);
ybaseTest = yBases(:,testFrame);
xRange = [minX maxX];
yRange = [minY maxY];

vMags = sqrt(vxTest.^2 + vyTest.^2);
vxTest(vMags>(mean(vMags) + 10*std(vMags)))=[];
vyTest(vMags>(mean(vMags) + 10*std(vMags)))=[];
xbaseTest(vMags>(mean(vMags) + 10*std(vMags)))=[];
ybaseTest(vMags>(mean(vMags) + 10*std(vMags)))=[];

vTest = binVectors(xbaseTest,ybaseTest,vxTest,vyTest,xRange,yRange,binSize, nThresh);

% Create vector matrix to be input into vectorFieldSparseInterpPatrick.m
% nx4 matrix of form [x0 y0 vx vy]n, (x0,y0) are base of vectors, (vx,vy) are the vectors themselves
% base is given by xCoords and yCoords, tips by xCoords+vx, similarly for y

adataVelFormatted = [xbaseTest ybaseTest vxTest vyTest];

% Make grid size mx2 of form [xg yg]m, grid size dr

dr = 1;
[xGrid, yGrid] = meshgrid(minX:dr:maxX, minY:dr:maxY);
gridMat = [xGrid(:) yGrid(:)];

threshold = 1;
d0 = 5*binSize; % Freedman et al 2016, S2
polygon = [];

% Interpolate test frame
%testVec = adataVelFormatted(:,testFrame);
testInterp = vectorFieldSparseInterpPatrick(adataVelFormatted, gridMat, threshold, d0, polygon);

interpedVX = reshape(testInterp(:,4), length(minX:dr:maxX), length(minX:dr:maxX));
interpedVY = reshape(testInterp(:,3), length(minY:dr:maxY), length(minY:dr:maxY));

% Now we compare the results of the interpolated field with the original one
% This part copied from plotsimdata.m

npoly=500;
adat=adata(:,:,testFrame);
subplot(2,2,1)
for k=1:1:npoly
    idx=find(adat(:,4)==k-1);
    X1=adat(idx(:),1);Y1=adat(idx(:),2);
    ind=rangesearch([X1 Y1],[X1(1) Y1(1)],10);
    h=plot(X1(ind{1}(:)),Y1(ind{1}(:)),'-r','MarkerSize',4,'MarkerFaceColor','r','LineWidth',2);
    hold on
end
str = sprintf('Original image, frame %d', testFrame);
xlim([-25, 25]), ylim([-25,25])
title(str)

subplot(2,2,4)
quiver(xGrid, yGrid, interpedVX, interpedVY)
title('Interpolated velocity field')
xlim([-25, 25]), ylim([-25,25])
subplot(2,2,3)
quiver(vTest(:,1), vTest(:,2), vTest(:,3), vTest(:,4))
title('Binned velocity field')
xlim([-25, 25]), ylim([-25,25])
subplot(2,2,2)
quiver(xbaseTest, ybaseTest, vxTest, vyTest)
title('Raw velocities')
xlim([-25, 25]), ylim([-25,25])

figure()
% Change NaN's to zeros first
interpedVX(isnan(interpedVX))=0;
interpedVY(isnan(interpedVY))=0;
div = divergence(xGrid, yGrid, interpedVX, interpedVY);
imagesc([-25,25],[-25,25],div')
colorbar
hold on
quiver(vTest(:,1), vTest(:,2), vTest(:,3), vTest(:,4),'k')
title(str)

toc