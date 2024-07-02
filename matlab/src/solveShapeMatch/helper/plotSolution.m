function fig = plotSolution(G, Va, Fa, Vb, Fb, not_matched_Fa, not_matched_Fb, shape2, varargin)
%PLOTSOLUTION Plots the solution G which was found with the shape-matching
%             algorithm. 
% Input Arguments:
%       - G:    binary vector containing the matchings in between shape A
%               and shape B. Either of size 3 * nFa * nFb or of size 
%               3 * nFa * nFb + 2 * 3 (nEa * nFb + nEb * nFa) + nVa * nFb +
%               nVb * nFa
%       - Va:   point coordinates from triangle a
%       - Fa:   triangulation from triangle a
%       - Vb:   point coordinates from triangle b
%       - Fb:   triangulation from triangle b

% settings
vertexSize = 40;
edgeThickness = 7;

degenerate = true;
plotDegenerate = true;

G = G > 0.5;
nFa = size(Fa, 1);
nFb = size(Fb, 1);
Ea = getEdgesNotClosed(Fa,0);
Eb = getEdgesNotClosed(Fb,0);
nEa = size(Ea,1);
nEb = size(Eb,1);
nVa = size(Va, 1);
nVb = size(Vb, 1);

if size(varargin{end},1)>0
    corr_infos = varargin{end};
    newComb = corr_infos('combNDeg');
    newComb = newComb{1};
    nNonDeg = size(newComb,1);
    FaComboNDeg = newComb(:,1:3);
    FbComboNDeg = newComb(:,4:6);
    degInfo = corr_infos('combDegA');
    degInfo = degInfo{1};
    FaComboDegA = degInfo(:,1:3);
    FbComboDegA = degInfo(:,4:6);
    degInfo = corr_infos('combDegB');
    degInfo = degInfo{1};
    FaComboDegB = degInfo(:,1:3);
    FbComboDegB = degInfo(:,4:6);
else
    [FaComboNDeg, FbComboNDeg] = getCombinations(Fa, Fb);
    nNonDeg = size(FaComboNDeg,1);

    if degenerate
        [FaComboDegA, FbComboDegA] = getDegenerateCombinations(Fa, Fb, 0);
        [FbComboDegB, FaComboDegB] = getDegenerateCombinations(Fb, Fa, 0);
    end
end

% nDegenerate = nNonDeg + 6 * nEa * nFb + 6* nEb * nFa...
%     + nVa * nFb + nVb * nFa;


% compute face combinations



fig = figure;
hold on;
cmap = hsv(size(FaComboNDeg,1));

colorsNDegA = zeros(size(FaComboNDeg));
for i=1:size(FaComboNDeg)
    curr_F = FaComboNDeg(i,:);
    curr_color = mean(Va(curr_F,:));
    curr_color(curr_color<0) = 0;
    curr_color(curr_color > 1) = 1;
    colorsNDegA(i,:) = curr_color;
end
%colorsNDegA = cmap;

% non-degenerate cases
g1 = G(1:nNonDeg);

ax1 = subplot(1,2,1);
plotFaces(Fa(not_matched_Fa, :), Va, 0.5 * ones(size(Fa(not_matched_Fa, :),1),3),0);
%plotFaces(Fa(not_matched_Fa, :), Va, 0.5 * ones(size(Fa(not_matched_Fa, :),1),3),0);
%[not_matched_Fb] = get_non_matched(FY_old, VY_old, FbComboNDeg(g1, :), Vb);
plotFaces(FaComboNDeg(g1, :), Va, colorsNDegA(g1,:),1);
xlim([min(Va(:,1)),max(Va(:,1))])
ylim([min(Va(:,2)),max(Va(:,2))])

hold on
ax2 = subplot(1,2,2);
plotFaces(Fb(not_matched_Fb, :), Vb, 0.5 * ones(size(Fb(not_matched_Fb, :),1),3),0);
plotFaces(FbComboNDeg(g1, :), Vb, colorsNDegA(g1,:),1);
%xlim([min(Va(:,1)),max(Va(:,1))])
%ylim([min(Va(:,2))-1,max(Va(:,2))-1])
hold on

% degenerate cases
if degenerate    
% degenerate cases in A
numDegA = size(FaComboDegA,1);
g2 = G(nNonDeg + 1:numDegA + nNonDeg);

colorsDegA = zeros(size(FaComboDegA));
for i=1:size(FaComboDegA)
    curr_F = FaComboDegA(i,:);
    curr_color = mean(Va(curr_F,:));
    curr_color(curr_color<0) = 0;
    curr_color(curr_color > 1) = 1;
    colorsDegA(i,:) = curr_color;
end

if plotDegenerate
    subplot(1,2,1)
    plotDegenerates(FaComboDegA,  Va, colorsDegA(g2, :),...
        g2, nFa, nFb, vertexSize, edgeThickness);
end
subplot(1,2,2)
plotFaces(FbComboDegA(g2, :), Vb, colorsDegA(g2, :),1);
% [not_matched_Fb] = get_non_matched(not_matched_Fb, VY_old, FbComboDegA(g2, :), Vb);
% plotFaces(not_matched_Fb, VY_old, 0.5 * ones(size(not_matched_Fb,1),3),0);


% degenerate cases in B
g3 = G(numDegA + nNonDeg + 1:end);

colorsDegB = zeros(size(FaComboDegB));
for i=1:size(FaComboDegB)
    curr_F = FaComboDegB(i,:);
    curr_color = mean(Va(curr_F,:));
    curr_color(curr_color<0) = 0;
    curr_color(curr_color > 1) = 1;
    colorsDegB(i,:) = curr_color;
end

subplot(1,2,1)
plotFaces(FaComboDegB(g3, :), Va, colorsDegB(g3,:),1);
if plotDegenerate
    subplot(1,2,2)
    plotDegenerates(FbComboDegB,  Vb, colorsDegB(g3,:),...
        g3, nFa, nFb, vertexSize, edgeThickness);
end

end


subplot(1,2,2)
view(3)
axis vis3d
%colorbar
axis equal
%colormap(hsv(100))

subplot(1,2,1)
view(3)
axis vis3d
%colorbar
axis equal
%colormap(hsv(100))
linkaxes([ax1 ax2],'xyz')



end


function plotFaces(F, V, cmap, not_matched)
patch('Faces', F, 'Vertices', V, 'FaceVertexCData', cmap, ...
    'FaceColor','flat','EdgeColor','w');
if not_matched
patch('Faces', F, 'Vertices', V, 'FaceVertexCData', cmap, ...
    'FaceColor','flat','EdgeColor','w');    
end
end

function plotDegenerates(Fa, Va, cmap, g, nFa, nFb, verSize, edgThick)
% points 
nVa = size(Va, 1);
numP2TriMatch = size(Fa,2);
gPoints = g(1:numP2TriMatch);
fPoints = Fa(gPoints, 1);
x = Va(fPoints, 1);
y = Va(fPoints, 2);
z = Va(fPoints, 3);
for i = 1:length(x)
    plot3(x(i), y(i), z(i), '.', 'MarkerSize', verSize, 'Color', cmap(i, :));
end
  
cmap = cmap(nnz(gPoints)+1:end, :);
Fa = Fa(numP2TriMatch+1:end, :);

% edges
gEdges = g(numP2TriMatch+1:end);
fEdges = Fa(gEdges, :);
if length(fEdges) == 0
    return;
end
% some index magic to extract the end points of the edges
idx1 = fEdges(:,2)==fEdges(:,3);
idx2 = idx1 == 0;
idx3 = (idx1 + idx2) ~= 2;
fEdgesIdx = [idx1 idx2 idx3];
fEdges = fEdges';
fEdges = fEdges(fEdgesIdx');
fEdges = reshape(fEdges, 2, nnz(fEdgesIdx)/2)';

% fEdges is now a numMatchedEdgesx2 vector

for i = 1:size(fEdges,1)
    line = [ Va(fEdges(i, 1), :);
             Va(fEdges(i, 2), :)];
    plot3(line(:, 1), line(:, 2), line(:, 3), ...
          '-', 'LineWidth', edgThick, 'Color', cmap(i, :));
end
end
