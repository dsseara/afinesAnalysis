function tensor = strain(xpos, ypos, dx, dy)
% strain takes a set of vectors on a grid and returns the associated symmetric strain tensor
%
% vBinned = binVectors(xBases, yBases, vx, vy, xRange, yRange, binSize)
% 
% INPUTS    xpos   : An [M 1] matrix with the bases of the vectors in the x direction
%           ypos   : An [N 1] matrix with the bases of the vectors in the y direction
%           dx     : An [M N] matrix with x displacements for each grid point
%           dy     : An [M N] matrix with y displacements for each grid point
%
%           NOTE   : If xpos and ypos are be given as a grid, they are converted to 
%                    a vector
%
% OUTPUTS   tensor : An [M N 2 2] matrix. For example, the symmetric strain tensor at
%                    position (x,y) is given by tensor(x,y,:,:)

% Created by Daniel Seara at 2017/01/17 13:32
% https://github.com/dsseara

% Check dimensionality of xpos and ypos
if sum(size(xpos)~=1) && sum(size(ypos)~=1)
    xpos = xpos(1,:);
    ypos = ypos(:,1);
elseif sum(size(xpos)~=1) || sum(size(ypos)~=1)
    error('Gave one grid and one vector for positions')
end

% Get spacing between grid points
hx = unique(diff(xpos));
hy = unique(diff(ypos));

A = [1 0 0 0; 0 0.5 0.5 0; 0 0 0 1]; %transformation matrix to go from derivs to symmetric
tensor = zeroes([size(dx) 2 2]);
[dxdx, dydx] = gradient(dx, hx);
[dxdy, dydy] = gradient(dy, hy);
for i=1:size(dx,1)
    for j=1:size(dx,2)
        E = (A*[dxdx(i,j); dydx(i,j); dxdy(i,j); dydy(i,j)])';
        e = [E(1) E(2); E(2) E(3)];
        tensor(i,j,:,:) = [xpos(i) ypos(j) e];
    end
end

end