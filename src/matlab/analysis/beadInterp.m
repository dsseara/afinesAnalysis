function nGrid = beadInterp(points,gridx,gridy)
% This function interpolates scattered points to the nearest point on a given grid,
% Assumes that only points within a circle with radius equal to half the diagonal between
% grid points are interpolated to each point
%
% INPUTS     points : An Nx2 array with bead positions (x,y)
%            gridx  : x coordinates of grid, the first output from meshgrid
%            gridy  : y coordinates of grid, the second output from meshgrid
%
% OUTPUTS    nGrid  : A grid of size(gridx) where each entry has the number of
%                     points interpolated to the corresponding grid point
%
% Created by Daniel Seara at 2017/02/28 17:28

nGrid = zeros(size(gridx));
N = size(points,1);
dx = unique(diff(gridx,1,2)); % Find grid spacing in x direction
dy = unique(diff(gridy,1,1)); % Find grid spacing in y direction

if numel(dx)~=1
    error('Found more than 1 grid size in x direction')
elseif numel(dy)~=1
    error('Found more than 1 grid size in y direction')
end

% Interpolation radius is half of the diagonal between grid points
radius = sqrt(dx.^2 + dy.^2)./2;

for ii=1:N
    point = points(ii,:);
    
    % Get a distance matrix
    distMat = sqrt((gridx-point(1)).^2 + (gridy-point(2)).^2);
    
    % Set all the values of distMat that are > r equal to 0
    distMat(distMat>radius) = 0;

    % Find the nonzero values, set them equal to 1/(# of nonzeros)
    % This lets 1/2 a particle be given to two points if it falls inside two interp regions
    distMat(distMat~=0) = 1/sum(sum(distMat~=0));
    
    nGrid = nGrid + distMat;
end