function projNonDegenerate = getNonDegenerateProjection(FX, FY)
% The product space for non degenerate cases is defined as follows:
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
%
%
% Consequently, the projection matrices have to be:
%
% projX \in IZ^{ numFacesX x numNonDegenerateCombos }
% projX = [  1  1  0  0  1  1  0  0  1  1  0  0       <- everywhere where
%                                                       face 1 2 3 is in
%                                                       Gamma goes a one
%            0  0  1  1  0  0  1  1  0  0  1  1]      <- same but for 4 5 6
% 
%       = repmat(repelem(speye(numFacesX), 1, numFacesY), 1, 3)
%
% 
% projY \in IZ^{ numFacesY x numNonDegenerateCombos }
% projY = [  1  0  1  0  1  0  1  0  1  0  1  0       <- for  7  8  9 
%            0  1  0  1  0  1  0  1  0  1  0  1]      <- for 10 11 12
%
%       = repmat(speye(numFacesY), 1, 3 * numFacesX)


numFacesX = size(FX, 1);
numFacesY = size(FY, 1);

projX = repmat(repelem(speye(numFacesX), 1, numFacesY), 1, 3);
projY = repmat(speye(numFacesY), 1, 3 * numFacesX);

projNonDegenerate = [projX; projY];

end
