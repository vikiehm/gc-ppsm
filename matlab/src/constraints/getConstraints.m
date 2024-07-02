function [constraintsMatrix, constraintsVector, product_E] = getConstraints(Va, Fa, Vb, Fb, degenerate, closed, corr_infos)
%getConstraints computes the constraints matrix and the constraints vector 
%               for a shape matching problem
%
% The constraints are defined as: 
%       A * Gamma = b
%
% where
% A = [                      del;          projX;          projY]
% b = [zeros(numProductEdges, 1); ones(numFacesX + numFacesY, 1)]
%
%
% Note: altDel option uses a different version of the del computation.
%       Read more on that in getDelAlt.

oldDel = false;
% this is needed in case we use mex functions
Fa = int32(Fa);
Fb = int32(Fb);
% compute the combinations
[FaComboNDeg, FbComboNDeg] = getCombinations(Fa, Fb);

if degenerate
    [FaComboDegA, FbComboDegA] = getDegenerateCombinations(Fa, Fb, closed);
    [FbComboDegB, FaComboDegB] = getDegenerateCombinations(Fb, Fa, closed);
    FaCombo = [FaComboNDeg; FaComboDegA; FaComboDegB];
    FbCombo = [FbComboNDeg; FbComboDegA; FbComboDegB];
else
    FaCombo = FaComboNDeg;
    FbCombo = FbComboNDeg;
end


%% del
[del, product_E] = getDel(FaCombo, FbCombo, Va, Fa, Vb, Fb, degenerate, closed, oldDel);

%% projection
[proj, projDegenerate] = getProjection(Fa, Fb, Va, Vb, degenerate);
projAll = [proj, projDegenerate];

constraintsMatrix = [del; projAll];

%% constraint vector
constraintsVector = getConstraintsVector(Va, Fa, Vb, Fb, degenerate, closed, oldDel);
end
