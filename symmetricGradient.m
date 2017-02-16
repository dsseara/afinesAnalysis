function [tensor,evecs,evals] = symmetricGradient(u,v,hx,hy)
% symmetricGradient takes a vector field on a grid and returns the symmetric part of the
% vector field's gradient
%
% [tensor, evecs, evals] = strain(u,v,hx,hy)
% 
% INPUTS         u : An [M N] matrix with x displacements for each grid point
%                v : An [M N] matrix with y displacements for each grid point
%               hx : grid spacing in x direction, default 1
%               hy : grid spacing in y direction, default 1
%           
%
% OUTPUTS   tensor : An [M N 2 2] matrix. For example, the symmetric strain tensor at
%                          position (x,y) is given by tensor(x,y,:,:)
%           evecs  : An [M N 2 2] matrix. For example, the columns of the matrix at
%                          evecs(x,y,:,:) are the eigenvectors of the symmetric tensor
%                          at tensor(x,y,:,:)
%           evals  : An [M N 2 2] matrix. For example, the diagaonal of the matrix at
%                          evals(x,y,:,:) are the eigenvalues of the symmetric tensor
%                          at tensor(x,y,:,:)
%
% Created by Daniel Seara at 2017/02/15 19:26
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

tensor = zeros([size(u) 2 2]);
evals  = zeros([size(u) 2 2]);
evecs  = zeros([size(u) 2 2]);
[dudx, dudy] = gradient(u, hx, hy);
[dvdx, dvdy] = gradient(v, hx, hy);

for i=1:size(u,1)
    for j=1:size(u,2)
        gradV = [dudx(i,j), dudy(i,j); dvdx(i,j), dvdy(i,j)];
        symTensor = (0.5) * (gradV + gradV');
        tensor(i,j,:,:) = symTensor;
        [evecs(i,j,:,:), evals(i,j,:,:)] = eig(symTensor);
    end
end

end