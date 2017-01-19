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
numTimeSteps = 2;

xyDisplacement = diff(adata(:,1:2, 1:numTimeSteps:end),1,3);

[numBeads, dim, numFrames] = size(xyDisplacement);
dt = numTimeSteps*uniquetol(diff(timestep));

xCoords = adata(:,1,1:numTimeSteps:end);
yCoords = adata(:,2,1:numTimeSteps:end);
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
binSize=2;
nThresh = 2;

% Make grid size mx2 of form [xg yg]m, grid size dr
dr = 2;
[xGrid, yGrid] = meshgrid(xRange(1):dr:xRange(2), yRange(1):dr:yRange(2));
gridMat = [xGrid(:) yGrid(:)];


interpRadius = 1; % radius of interpolation region
d0 = 5*binSize; % Freedman et al 2016, S2, width of Gaussian to weight
polygon = [];

% Initialize arrays for storage of strain and Strain Rate time series
strainSeries     = zeros(numFrames,1);
strainRateSeries = zeros(numFrames,1);

adataVelInterped = NaN(size(xGrid) dim numFrames);

parfor i=1:numFrames
    disp(i)
    vxFrame = vx(:,i);
    vyFrame = vy(:,i);
    xbaseFrame = xBases(:,i);
    ybaseFrame = yBases(:,i);
    
    vBinned  = binVectors(xbaseFrame,ybaseFrame,vxFrame,vyFrame,xRange,yRange,binSize,nThresh);
    %drBinned = binVectors(xbaseFrame,ybaseFrame,dxFrame,dyFrame,xRange,yRange,binSize,nThresh);

    % Interpolate frame
    vInterp = vectorFieldSparseInterpPatrick(vBinned,  gridMat, interpRadius, d0, polygon);
    %rInterp = vectorFieldSparseInterpPatrick(drBinned, gridMat, interpRadius, d0, polygon);

    interpedVX = reshape(vInterp(:,3), length(xGrid(1,:)), length(xGrid(:,1)));
    interpedVY = reshape(vInterp(:,4), length(yGrid(1,:)), length(yGrid(:,1)));
    interpedVX(isnan(interpedVX)) = 0;
    interpedVY(isnan(interpedVY)) = 0;

    adataVelInterped(:,:,1,i) = interpedVX;
    adataVelInterped(:,:,2,i) = interpedVY;
    save([pwd,'/interpedData.mat'])
    % npoly=500;
    % adat=adata(:,:,actualFrames(i));
    % for k=1:1:npoly
    %     idx=find(adat(:,4)==k-1);
    %     X1=adat(idx(:),1);Y1=adat(idx(:),2);
    %     ind=rangesearch([X1 Y1],[X1(1) Y1(1)],10);
    %     h=plot(X1(ind{1}(:)),Y1(ind{1}(:)),'-r','MarkerSize',4,'MarkerFaceColor','r','LineWidth',2);
    %     hold on
    % end
    % quiver(xGrid, yGrid, interpedVX, interpedVY,'k')
    % str = sprintf('Interpolated overlay, time step %0.2f', i*dt);
    % xlim([-25, 25]), ylim([-25,25])
    % title(str)
    % %saveas(gcf,['Overlay/',num2str(i)],'tif')
    % hold off
    % pause
    % close all


    % interpedDX = reshape(rInterp(:,3), length(xGrid(1,:)), length(xGrid(:,1)));
    % interpedDY = reshape(rInterp(:,4), length(yGrid(1,:)), length(yGrid(:,1)));
    % interpedDX(isnan(interpedDX)) = 0;
    % interpedDY(isnan(interpedDY)) = 0;

%    [strainTensor, strainEvec, strainEvals] = strain(interpedDX, interpedDY, dr);
    [strainRateTensor, strainRateEvecs, strainRateEvals] = strain(interpedVX, interpedVY, dr);
%    strainSeries(i) = sum(sum(sum(sum(strainEvals))));
    strainRateSeries(i) = mean(mean(sum(sum(strainRateEvals,3),4)));
end % end loop over all the frames


% Now we compare the results of the interpolated field with the original one
% This part copied from plotsimdata.m

% npoly=500;
% adat=adata(:,:,testFrame);
% subplot(2,2,1)
% for k=1:1:npoly
%     idx=find(adat(:,4)==k-1);
%     X1=adat(idx(:),1);Y1=adat(idx(:),2);
%     ind=rangesearch([X1 Y1],[X1(1) Y1(1)],10);
%     h=plot(X1(ind{1}(:)),Y1(ind{1}(:)),'-r','MarkerSize',4,'MarkerFaceColor','r','LineWidth',2);
%     hold on
% end
% str = sprintf('Original image, frame %d', testFrame);
% xlim([-25, 25]), ylim([-25,25])
% title(str)

% subplot(2,2,4)
% quiver(xGrid, yGrid, interpedVX, interpedVY)
% title('Interpolated velocity field')
% xlim([-25, 25]), ylim([-25,25])
% subplot(2,2,3)
% quiver(vTest(:,1), vTest(:,2), vTest(:,3), vTest(:,4))
% title('Binned velocity field')
% xlim([-25, 25]), ylim([-25,25])
% subplot(2,2,2)
% quiver(xbaseTest, ybaseTest, vxTest, vyTest)
% title('Raw velocities')
% xlim([-25, 25]), ylim([-25,25])

% figure()
% % Change NaN's to zeros first
% interpedVX(isnan(interpedVX))=0;
% interpedVY(isnan(interpedVY))=0;
% div = divergence(xGrid, yGrid, interpedVX, interpedVY);
% pcolor(xGrid, yGrid, div)
% shading interp;
% axis image;
% colorbar
% hold on
% quiver(xGrid, yGrid, interpedVX, interpedVY,'k')
% title(str)

toc