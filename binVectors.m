function vBinned = binVectors(xBase, yBase, vx, vy, xRange, yRange, binSize, threshold)
% binVectors takes a vector field and bins the vectors according to their starting position.
% For each bin, the average vector and starting position is returned. 
% First removes any vectors that are larger than mean(magnitudes) + 10*std(magnitudes)
% 
% vBinned = binVectors(xBases, yBases, vx, vy, xRange, yRange, binSize)
% 
% INPUTS    xBase     : (Nx1) vector of x-coordinate of base of vectors                  |
%           yBase     : (Nx1) vector of y-coordinate of base of vectors                  | These four are identical
%           vx        : (Nx1) vector of x-component of vector at corresponding (x,y)     | to quiver input
%           vy        : (Nx1) vector of y-component of vector at corresponding (x,y)     |
%           xRange    : size of velocity field area, a (1x2) vector of form [xMin xMax]
%           yRange    : size of velocity field area, a (1x2) vector of form [xMin xMax]
%           binSize   : size of bin, a scalar (i.e. only square bins for now), same units
%                       as xRange and yRange
%           threshold : if number of vectors in bin < threshold, count as zero vectors in that bin
% 
% OUTPUTS   vBinned : (Mx4) matrix of binned vector of the form [x0' y0' vx' vy'],
%                     similar format to V.
%                     M = (xMax - xMin)/dx * (yMax - yMin)/dy
%
% Created by Daniel Seara at 2017/01/11 16:08
% https://github.com/dsseara


vBinned = [];

% Remove huge vectors
vMags = sqrt(vx.^2 + vy.^2);
vx(vMags>(mean(vMags) + 10*std(vMags)))=[];
vy(vMags>(mean(vMags) + 10*std(vMags)))=[];
xBase(vMags>(mean(vMags) + 10*std(vMags)))=[];
yBase(vMags>(mean(vMags) + 10*std(vMags)))=[];

binedgesX = xRange(1):binSize:xRange(2);
binedgesY = yRange(1):binSize:yRange(2);

for i = 1:numel(binedgesX)-1
    xEdge = binedgesX(i);
    indX = xBase >= xEdge & xBase < xEdge+binSize;
    for j = 1:numel(binedgesY)-1
        yEdge = binedgesY(j);
        indY = yBase >= yEdge & yBase < yEdge+binSize;
        inThisBin = indX & indY;
        
        if sum(inThisBin)<threshold % do nothing if there's not enough vectors in the bin
        else
            vBinned = [vBinned; mean(xBase(inThisBin)) mean(yBase(inThisBin)) mean(vx(inThisBin)) mean(vy(inThisBin))];
        end
    end
end

if isempty(vBinned)
    vBinned = zeros(1,4);
end

%binnedVecs = [binnedVecs; repmat(xEdge, 1,1,numFrames) repmat(yEdge, 1,1,numFrames) reshape(sum(ind),1,1,numFrames)];