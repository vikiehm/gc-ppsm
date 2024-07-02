function [FaCombo, FbCombo] = getCombinations(Fa, Fb)
% Consider the face matrix of shape A  and the one of shape B
%      Fa = [1 2 3;         Fb = [ 7  8  9; 
%            4 5 6];              10 11 12];
% Than the matchings are definied as: 
%       Matchings = [[1 2 3] [ 7  8  9]]
%                   [[1 2 3] [10 11 12]]
%                   [[4 5 6] [ 7  8  9]]
%                   [[4 5 6] [10 11 12]]
%                   [[2 3 1] [ 7  8  9]]
%                   [[2 3 1] [10 11 12]]
%                   [[5 6 4] [ 7  8  9]]
%                   [[5 6 4] [10 11 12]]
%                   [[3 1 2] [ 7  8  9]]
%                   [[3 1 2] [10 11 12]]
%                   [[6 4 5] [ 7  8  9]]
%                   [[6 4 5] [10 11 12]]
% returns column vector for indexing each of Fa and Fb

rot = [2 3 1];

nFa = size(Fa, 1);
nFb = size(Fb, 1);

% combine all faces in A with all faces in B
faRowIdx = repelem(1:nFa, nFb);
fbRowIdx = repmat(1:nFb, 1, nFa)';

% first: with no rotation in A
FaCombo = Fa(faRowIdx, :);
FbCombo = Fb(fbRowIdx, :);

% second: we rotate A once
Fa = Fa(:, rot);
FaCombo = [FaCombo; Fa(faRowIdx, :)];
FbComboTemp = FbCombo;
FbCombo = [FbCombo; FbCombo];

% third: we rotate A twice
Fa = Fa(:, rot);
FaCombo = [FaCombo; Fa(faRowIdx, :)];
FbCombo = [FbCombo; FbComboTemp];

end