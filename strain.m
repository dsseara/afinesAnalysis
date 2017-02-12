function [strainTensor,strainEvecs,strainEvals] = strain(dx,dy,hx,hy)
% strain takes a set of vectors on a grid and returns the associated symmetric strain tensor
%
% [strainTensor, strainEvecs, strainEvals] = strain(dx,dy,hx,hy)
% 
% INPUTS              dx : An [M N] matrix with x displacements for each grid point
%                     dy : An [M N] matrix with y displacements for each grid point
%                     hx : grid spacing in x direction, default 1
%                     hy : grid spacing in y direction, default 1
%           
%
% OUTPUTS   strainTensor : An [M N 2 2] matrix. For example, the symmetric strain tensor at
%                          position (x,y) is given by tensor(x,y,:,:)
%           strainEvecs  : An [M N 2 2] matrix. For example, the columns of the matrix at
%                          strainEvecs(x,y,:,:) are the eigenvectors of the symmetric strain tensor
%                          at tensor(x,y,:,:)
%           strainEvals  : An [M N 2 2] matrix. For example, the diagaonal of the matrix at
%                          strainEVals(x,y,:,:) are the eigenvalues of the symmetric strain tensor
%                          at tensor(x,y,:,:)
%
% Created by Daniel Seara at 2017/01/17 13:32
% https://github.com/dsseara


if nargin<2
    error('Need two scalar fields to describe a vector field \n')
elseif nargin==2
    hx=1;
    hy=1;
elseif nargin==3
    hy=hx;
elseif nargin>4
    error('Too many inputs \n');
end

strainTensor = zeros([size(dx) 2 2]);
strainEvals  = zeros([size(dx) 2 2]);
strainEvecs  = zeros([size(dx) 2 2]);
[dxdx, dxdy] = gradient(dx, hx, hy);
[dydx, dydy] = gradient(dy, hx, hy);

for i=1:size(dx,1)
    for j=1:size(dx,2)
        gradV = [dxdx(i,j), dxdy(i,j); dydx(i,j), dydy(i,j)];
        symStrain = (0.5) * (gradV + gradV');
        strainTensor(i,j,:,:) = symStrain;
        [strainEvecs(i,j,:,:), strainEvals(i,j,:,:)] = eig(symStrain);
    end
end

end