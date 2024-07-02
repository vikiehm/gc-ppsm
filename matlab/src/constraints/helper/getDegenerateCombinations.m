function [FaCombo, FbCombo] = getDegenerateCombinations(Fa, Fb, closed)
% Consider the face matrix of shape A  and the one of shape B
%      Fa = [1 2 3;         Fb = [ 7  8  9; 
%            4 5 6];              10 11 12];
%
% => EdgesA = [1 2;
%              2 3;
%              3 1;
%              4 5;
%              5 6;
%              6 4];
% => From one edge we can generate two different degenerate triangles: 
% e.g edge = [1 2] 
%     tri1 = [1 1 2]
%     tri2 = [1 2 2]
%
%
%       Matchings = 
% The degenerate matchings in A are definied as:
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
% ...
% Than the degenerate matchings in B are definied as:
%   ... analogously to A
%
%
% returns column vector for indexing each of Fa and Fb

nFa = size(Fa, 1);
nFb = size(Fb, 1);
% only if shape A is closed
if closed
    nEa = nFa * 3/2; 
    assert(floor(nEa) == nEa, 'Shape Closed?');
    nVa = max(max(Fa)); 

else
   [Ea, ~] = getEdgesNotClosed(Fa, closed);
    %Ea = getEdges(Fa, closed);
    nEa = size(Ea,1);
    nVa = size(unique(Fa),1);
end

nDegenerateACombos = nVa * nFb + 2 * nEa * nFb * 3;

% FaCombo
FaCombo = zeros(nDegenerateACombos, 3);
%FbCombo = zeros(nFb * nFa * 42, 3);

% triangle to vertex
FaCombo(1:nVa * nFb, :) = repmat(repelem((1:nVa)', nFb, 1), 1, 3);  

% triangle to edge
FaCombo(nVa * nFb + 1: end, :) = getTriangle2EdgeMatching(Fa, nFb, closed);

% FbCombo
rot1 = [2 3 1];
rot2 = [3 1 2];
FbCombo = [repmat(Fb, nVa, 1);...              % triangle to vertex
           repmat(Fb, nEa * 2, 1);...          % triangle to edge 1st rot
           repmat(Fb(:, rot1), nEa * 2, 1);... % triangle to edge 2nd rot
           repmat(Fb(:, rot2), nEa * 2, 1)];   % triangle to edge 3rd rot
end


function FaCombo = getTriangle2EdgeMatching(Fa, nFb, closed)
Ea = getEdgesNotClosed(Fa, closed);

edge2Tri1 = [1 1 2];
edge2Tri2 = [1 2 2];

FaCombo = [Ea(:, edge2Tri1);
           Ea(:, edge2Tri2);];
       
FaCombo = repelem(FaCombo, nFb, 1);
FaCombo = repmat(FaCombo, 3, 1);
end