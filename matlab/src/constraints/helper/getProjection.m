function [proj, projDegenerate]  = getProjection(FX, FY, VX, VY, degenerate)
%getProjection computes the projection part of the constraints

proj = getNonDegenerateProjection(FX, FY);
projDegenerate = [];
if degenerate 
    [projX_XDeg, projY_XDeg] = getDegenerateProjection(FX, FY, VX);
    [projY_YDeg, projX_YDeg] = getDegenerateProjection(FY, FX, VY);
     projX = [projX_XDeg, projX_YDeg];
     projY =  [projY_XDeg, projY_YDeg];

     projDegenerate = [projX; projY];            
end
end


function [projX, projY] = getDegenerateProjection(FX, FY, VX)
% The product space for degenerate cases in X is defined as follows:
%      Fx = [1 2 3;         Fy = [ 7  8  9; 
%            4 5 6];              10 11 12];
%
%       Matchings = 
% The degenerate matchings in X are definied as:
% Triangle to vertex
%                   [[1 1 1] [ 7  8  9]]
%                   [[1 1 1] [10 11 12]]
%                   [[2 2 2] [ 7  8  9]]
%                   [[2 2 2] [10 11 12]]
%                   [[3 3 3] [ 7  8  9]]
%                   [[3 3 3] [10 11 12]]
%                   [[4 4 4] [ 7  8  9]]
%                   [[4 4 4] [10 11 12]]
%                   [[5 5 5] [ 7  8  9]]
%                   [[5 5 5] [10 11 12]]
%                   [[6 6 6] [ 7  8  9]]
%                   [[6 6 6] [10 11 12]]
% Triangle to edge (X)
% 1. Rotation
%                   [[1 1 2] [ 7  8  9]]
%                   [[1 1 2] [10 11 12]]
%                   [[2 2 3] [ 7  8  9]]
%                   [[2 2 3] [10 11 12]]
%                   [[3 3 1] [ 7  8  9]]
%                   [[3 3 1] [10 11 12]]
%                   [[4 4 5] [ 7  8  9]]
%                   [[4 4 5] [10 11 12]]
%                   [[5 5 6] [ 7  8  9]]
%                   [[5 5 6] [10 11 12]]
%                   [[6 6 4] [ 7  8  9]]
%                   [[6 6 4] [10 11 12]]
%
%                   [[1 2 2] [ 7  8  9]]
%                   [[1 2 2] [10 11 12]]
%                   [[2 3 3] [ 7  8  9]]
%                   [[2 3 3] [10 11 12]]
%                   [[3 1 1] [ 7  8  9]]
%                   [[3 1 1] [10 11 12]]
%                   [[4 5 5] [ 7  8  9]]
%                   [[4 5 5] [10 11 12]]
%                   [[5 6 6] [ 7  8  9]]
%                   [[5 6 6] [10 11 12]]
%                   [[6 4 4] [ 7  8  9]]
%                   [[6 4 4] [10 11 12]]
% 2. Rotation
%                   [[1 1 2] [ 8  9  7]]
%                   [[1 1 2] [11 12 10]]
%                   [[2 2 3] [ 8  9  7]]
% ...
% 3. Rotation
%                   [[1 1 2] [ 9  7  8]]
%                   [[1 1 2] [12 10 11]]
%                   [[2 2 3] [ 9  7  8]]
%
% The product space for non degenerate cases in X is defined as follows:
% * All triangles of X are degenerate => all projections of X are zero
%   projX \in IZ^{ numFacesX x numDegenerateCombos }
%   projX = [zeros(numFacesX, numTriangle2Vertex),
%                zeros(numFacesX, 3 * numTriangleToEdge)]
% * The triangles in Y are just repeated several times:
%   projY \in IZ^{ numFacesX x numDegenerateCombos }
%   projY = [repmat(speye(numFacesY), 1, 3 * numFacesX), 
%            repmat(speye(numFacesY), 1,  18 * numFacesX)]

numFacesX = length(FX);
numFacesY = length(FY);
numVerticesX = size(VX,1);
numTriangle2Vertex  = numVerticesX * numFacesY;
[EX,~] = getEdgesNotClosed(FX, 0);
numTriangleToEdge   = 6 * numFacesY * size(EX,1);


projX = [sparse(numFacesX, numTriangle2Vertex), ...
         sparse(numFacesX, numTriangleToEdge)];

projY = [repmat(speye(numFacesY), 1, numVerticesX), ...
         repmat(speye(numFacesY), 1,  6*size(EX,1))];



end
