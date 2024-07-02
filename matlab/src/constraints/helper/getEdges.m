function [E, edgeLocationsInF] = getEdges(F, closed)
%getEdges returns a list with the edges contained in the face matrix F
%  closed: if shape is closed
% * If two output arguments are specified additionally the position of the
%   edge and its opposite directed edged are stored in a numEdges x 2 
%   matrix.
% * Algorithm: It traverses the matrix F and extracts the edges contained 
%   in F. For closed shapes each edge appears two times in F. => the 
%   function just adds one orientation.
% 
%
% Note: The shape has to be closed for this algorithm to work!

numFaces = size(F, 1);
numEdges = numFaces * 3 / 2;
if closed
    assert(floor(numEdges) == numEdges, 'Is your shape closed?');
else
    % remove to much edges in the end
    numEdges = numFaces * 3;
end
% prealloc
E = zeros(numEdges, 2, 'int32');
edgeLocationsInF = -ones(numEdges, 2, 'int32');

% add the edges of the first triangle manually, otherwise the loop
% doesn't work
E(1, :) = F(1, [1 2]);
edgeLocationsInF(1, 1) = 1;
E(2, :) = F(1, [2 3]);
edgeLocationsInF(2, 1) = 1;
E(3, :) = F(1, [3 1]);
edgeLocationsInF(3, 1) = 1;
numEdgesAdded = 3;

%% loop through all other triangles
for f = 2:numFaces
    edge1 = F(f, [1 2]); found1 = false;
    edge2 = F(f, [2 3]); found2 = false;
    edge3 = F(f, [3 1]); found3 = false;
    inv = [2 1];
    minusEdge1 = edge1(inv);
    minusEdge2 = edge2(inv);
    minusEdge3 = edge3(inv);
    
    % loop through all already added edges
    for e = 1:numEdgesAdded
        
        % edge 1
        if all(E(e, :) == edge1)
            found1 = true; 
        end
        if all(E(e, :) == minusEdge1)
            edgeLocationsInF(e, 2) = f;
            found1 = true; 
        end
        
        % edge 2
        if all(E(e, :) == edge2) 
            found2 = true; 
        end
        if all(E(e, :) == minusEdge2)
            edgeLocationsInF(e, 2) = f;
            found2 = true; 
        end
        
        % edge 3
        if all(E(e, :) == edge3)
            found3 = true; 
        end
        if all(E(e, :) == minusEdge3)
            edgeLocationsInF(e, 2) = f;
            found3 = true; 
        end
        
    end
    
    
    %% Add edges to E matrix if new edges were found
    if numEdgesAdded > numEdges
        % this should never happen
        error('Is your shape closed?')
    end
    
    if ~found1
        numEdgesAdded = numEdgesAdded + 1;
        E(numEdgesAdded, :) = edge1;
        edgeLocationsInF(numEdgesAdded, 1) = f;
    end

    if ~found2
        numEdgesAdded = numEdgesAdded + 1;
        E(numEdgesAdded, :) = edge2;
        edgeLocationsInF(numEdgesAdded, 1) = f;
    end

    if ~found3
        numEdgesAdded = numEdgesAdded + 1;
        E(numEdgesAdded, :) = edge3;
        edgeLocationsInF(numEdgesAdded, 1) = f;
    end
    
    
end

if numEdgesAdded ~= numEdges
    if closed
        error("Extracted less edges than expected. Is your mesh closed?");
    else
        E = E(1:numEdgesAdded,:);
        edgeLocationsInF = edgeLocationsInF(1:numEdgesAdded,:);
    end
end

% mesh is only watertight if each does appear exactly two times (in
% opposite directions)
if closed
    if any(edgeLocationsInF == -1)
        error("Mesh is not closed!");
    end
end

end