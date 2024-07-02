function [Del,E] = getDel(FXCombo, FYCombo, VX, FX, VY, FY, degenerate, closed, old)
%getDel returns the del part of the constraints matrix
%
% We can construct the del part as follows:
% iterate through all triangles in the product space => column pos in del
%   iterate through all edges in the in the edge product space
%       if oriented as in edge list => 1 in del matrix
%       if not oriented as in edge list => -1 in del matrix
% Since this is very computational demanding ( O(numProductFaces^2) ), we 
% propose an approach that cleverly utilizes the structure of the product 
% space and the location of the edges of each shape within their face 
% matrices.
%
% The del matrix can be seen as follows:
% 
%                                 numProductFaces
%         ┌──────────────────────────────────────────────────────────────┐
%         │                                                              │
%         │                                                              │
%  num    │                                                              │
% Product │                                                              │
%  Edges  │                                                              │
%         │                                                              │
%         │                                                              │
%         └──────────────────────────────────────────────────────────────┘
% We now can split this into different sections:
%
%     non degenerate faces
%              | vertex X to triangle Y         vertex Y to triangleX
%              |      |    edgesX to triangleY     |   edgesY to triangleX
%              |      |             |              |           |
%          ┌──────|──────|─────────────────────|──────|──────────────────┐
% non deg  │      |      |                     |      |                  │
%   Edges  │      |      |                     |      |                  │
% ---------- ----------------------------------------------------------- -
% deg Edges│      |      |                     |      |                  │
%  of X    │      |      |                     |      |                  │
% ---------- ----------------------------------------------------------- -
% deg Edges│      |      |                     |      |                  │
%  of      │      |      |                     |      |                  │
%          └─────────────────────────────────────────────────────────────┘
%
% We know: 
% * non deg Edges can't appear in vertex to triangle product faces
% * deg edges in Z can only appear in vertex Z to triangle and in edgesZ to
%   triangle product faces
% * the ordering of the faces in the product space is according to their
%   order in the respective face matrices 
% * the ordering of the edges in the edges product space is dependent on
%   the order of the faces (if no sorting is performed)
%
% => combining these facts and cleverly indexing only the relevant faces by
% combing the information of the location of the edges we can reduce the
% complexity of the getDel function to O(numProductFaces)
%
%
% REFERENCES:
% (1) WINDHEUSER, Thomas, et al. 
% Large‐scale integer linear programming for orientation preserving 3d 
% shape matching. 
% In: Computer Graphics Forum. Oxford, UK: Blackwell Publishing Ltd, 
% 2011. S. 1471-1480.

% checkEdge = @checkEdgeNew;
% if old
%      checkEdge = @checkEdgeOld;
% end

% get edges and the corresponding locations
[edgesX, locEdgesXinFX] = getEdgesNotClosed(FX, closed);
[edgesY, locEdgesYinFY] = getEdgesNotClosed(FY, closed);

onlyComputeProductEdges = false;
numFaces = size(FXCombo, 1)
numVerticesX = size(VX, 1);
numVerticesY = size(VY, 1);
numFacesX = size(FX, 1);
numFacesY = size(FY, 1);
numFacesXxNumFacesY = numFacesX * numFacesY;
if closed
    numEdgesX = numFacesX * 3 / 2;
    numEdgesY = numFacesY * 3 / 2;
else
    numEdgesX = size(edgesX,1);
    numEdgesY = size(edgesY,1);
end
numEdgesNonDegenerate = numEdgesX * numEdgesY;
numEdgesYNonDegenerate = numEdgesNonDegenerate + numVerticesX * numEdgesY;
if degenerate 
    numEdges =  numEdgesNonDegenerate + numVerticesX * numEdgesY + ...
                numEdgesX * numVerticesY;
else
    numEdges = numEdgesNonDegenerate;
end

eOld = 0;
if old
    eOld = numEdgesNonDegenerate;
    numEdges = numEdges + numEdgesNonDegenerate;
end


%% Construct edge product space and edge product space to edge translators

E = [repelem(edgesX, numEdgesY, 1), repmat(edgesY, numEdgesX, 1)];
eToEdgesXTranslator = repelem((1:numEdgesX)', numEdgesY, 1);
eToEdgesYTranslator = repmat( (1:numEdgesY)', numEdgesX, 1);
old = false;
if old 
    disp("old")
    E = [E; repelem(edgesX(:, [2 1]), numEdgesY, 1), repmat(edgesY, numEdgesX, 1)];
    eToEdgesXTranslator = [eToEdgesXTranslator; repelem((1:numEdgesX)', numEdgesY, 1)];
    eToEdgesYTranslator = [eToEdgesYTranslator; repmat( (1:numEdgesY)', numEdgesX, 1)];
end

if degenerate
    disp("degenerate")
    degEdgesX = repmat((1:numVerticesX)', 1, 2);
    degEdgesY = repmat((1:numVerticesY)', 1, 2);
    degEdgs = [repelem(degEdgesX, numEdgesY, 1), repmat(edgesY, numVerticesX, 1);
               repelem(edgesX, numVerticesY, 1), repmat(degEdgesY, numEdgesX, 1)];
    E = [   E; 
            degEdgs];
    size(E)

    if onlyComputeProductEdges
        Del = sparse(size(E,1), size(FXCombo,1));
        return
    end

        
    eToEdgesXTranslator = [eToEdgesXTranslator; repelem((1:numVerticesX)', numEdgesY, 1)];
    eToEdgesXTranslator = [eToEdgesXTranslator; repelem((1:numEdgesX)', numVerticesY, 1)]; 
    eToEdgesYTranslator = [eToEdgesYTranslator; repmat((1:numEdgesY)', numVerticesX, 1)];
    eToEdgesYTranslator = [eToEdgesYTranslator; repmat((1:numVerticesY)', numEdgesX, 1)];
end

%
oldVertexIdx = -1;
perIdxes = zeros(50, 2);
perNumIdxes = zeros(1, 2);

% prealloc the sparse array for speed (in each column there should be 
% 3 nnz elements)
nnzDel = 3 * numFaces;

del = struct('i', zeros(nnzDel, 1), 'j', zeros(nnzDel, 1), 'val', zeros(nnzDel, 1), 'length', 0);
%del = spalloc(length(E), numFaces, nnzDel);

%% Iterate through all edges in E

% non-degenerate edges
for e = 1:numEdgesNonDegenerate + eOld
    searchInNonDegenerateFaces(e);
    if degenerate
        searchInEdgesX2TriangleY(e);
        searchInEdgesYToTriangleX(e);
    end
     
end

% degenerate edges
if degenerate
    
    % reset persistent variables in this function
    findVerticesInEdgesMatrix(0, 1, 0); 
    
    % degenerate edges in X
    for e = eOld + numEdgesNonDegenerate+1 : eOld + numEdgesYNonDegenerate
        searchInVertexX2TriangleY(e);
        searchDegEdgesXInEdgesX2TriangleY(e);
    end

    % reset persistent variables in this function
    findVerticesInEdgesMatrix(0, 1, 0);
    
    % degenerate edges in Y
    for e = eOld + numEdgesYNonDegenerate+1:numEdges
        searchInVertexY2TriangleX(e);
        searchDegEdgesYInEdgesY2TriangleX(e);
    end
end


if closed== 1 && del.length ~= nnzDel
    disp(del.length)
    error('not enough elements added')
end
Del = sparse(del.i(1:del.length), del.j(1:del.length), del.val(1:del.length));


%% Search functions for the different part of the del matrix

function searchInNonDegenerateFaces(e)
    % this i-loop accounts for the different rotations of the triangles
    for i = 0:3-1
        f = i * numFacesXxNumFacesY;
        % each product edge may be in one of 4 possible triangles
        
        f11 = f + (locEdgesXinFX(eToEdgesXTranslator(e), 1) - 1) * numFacesY...
                + locEdgesYinFY(eToEdgesYTranslator(e), 1);
        val = checkFace(e, f11);

        if val
            del.length = del.length + 1;
            del.i(del.length) = e;
            del.j(del.length) = f11;
            del.val(del.length) = val;
        end
        
        f12 = f + (locEdgesXinFX(eToEdgesXTranslator(e), 1) - 1) * numFacesY...
                + locEdgesYinFY(eToEdgesYTranslator(e), 2);
        val = checkFace(e, f12);
        if val
            del.length = del.length + 1;
            del.i(del.length) = e;
            del.j(del.length) = f12;
            del.val(del.length) = val;
        end
        
        f21 = f + (locEdgesXinFX(eToEdgesXTranslator(e), 2) - 1) * numFacesY...
                + locEdgesYinFY(eToEdgesYTranslator(e), 1);
        val = checkFace(e, f21);
        if val
            del.length = del.length + 1;
            del.i(del.length) = e;
            del.j(del.length) = f21;
            del.val(del.length) = val;
        end
        
        f22 = f + (locEdgesXinFX(eToEdgesXTranslator(e), 2) - 1) * numFacesY...
                + locEdgesYinFY(eToEdgesYTranslator(e), 2);
        val = checkFace(e, f22);
        if val
            del.length = del.length + 1;
            del.i(del.length) = e;
            del.j(del.length) = f22;
            del.val(del.length) = val;
        end
    end
    
end


function searchInVertexX2TriangleY(e)
    offset = 3 * numFacesXxNumFacesY;
    f = offset + (eToEdgesXTranslator(e)-1) * numFacesY;
    
    % first orientation of the edge
    f1 = f + locEdgesYinFY(eToEdgesYTranslator(e), 1);
    
    val = checkFace(e, f1);
    if val
        del.length = del.length + 1;
        del.i(del.length) = e;
        del.j(del.length) = f1;
        del.val(del.length) = val;
    end
    
    % second orientation of the edge
    f2 = f + locEdgesYinFY(eToEdgesYTranslator(e), 2);
    val = checkFace(e, f2);
    if val
        del.length = del.length + 1;
        del.i(del.length) = e;
        del.j(del.length) = f2;
        del.val(del.length) = val;
    end
    
end


function searchInEdgesX2TriangleY(e)
    offset = 3 * numFacesXxNumFacesY + numVerticesX * numFacesY;
    
    % the i loop accounts for the possible rotations
    for i = 0:3-1 
        rotPos = i * 2 * numEdgesX * numFacesY;
        
        % the ii loop accounts for different possibilities the triangles
        % can be constructet from one edge
        for ii = 0:2-1 
            edge2TriPos = ii * numEdgesX * numFacesY;
            
            f = offset + edge2TriPos + rotPos +...
                    (eToEdgesXTranslator(e) - 1) * numFacesY;
                
            % first orientation of the edge
            f1 = f + locEdgesYinFY(eToEdgesYTranslator(e), 1);
            val = checkFace(e, f1);
            if val
                del.length = del.length + 1;
                del.i(del.length) = e;
                del.j(del.length) = f1;
                del.val(del.length) = val;
            end
            
            % second orientation of the edge
            f2 = f + locEdgesYinFY(eToEdgesYTranslator(e), 2);
            val = checkFace(e, f2);
            if val
                del.length = del.length + 1;
                del.i(del.length) = e;
                del.j(del.length) = f2;
                del.val(del.length) = val;
            end
        end
    end
end


function searchDegEdgesXInEdgesX2TriangleY(e)
    offset = 3 * numFacesXxNumFacesY + numVerticesX * numFacesY;
    
    % we don't know how often a vertex appears in the corresponding edge
    % matrix of a shape => this functions searches for the vertices
    [idxes, numIdxes] = findVerticesInEdgesMatrix(eToEdgesXTranslator(e), 0, 0);
    
    % the i loop accounts for the possible rotations
    for i = 0:3-1 
        rotPos = i * 2 * numEdgesX * numFacesY;
        
        % the vertex index can either be in the first or the second element
        % of an edge => k == 1 first element; k == 2 second element
        for k = 1:2
            for idx = 1:numIdxes(k)
                edge2TriPos = (idxes(idx, k)-1) * numFacesY ...
                    + (k-1) * numEdgesX * numFacesY;

                % first orientation of the edge
                f = offset + edge2TriPos + rotPos;
                f1 = f + locEdgesYinFY(eToEdgesYTranslator(e), 1);
                val = checkFace(e, f1);
                if val
                    del.length = del.length + 1;
                    del.i(del.length) = e;
                    del.j(del.length) = f1;
                    del.val(del.length) = val;
                end
                
                % second orientation of the edge
                f2 = f + locEdgesYinFY(eToEdgesYTranslator(e), 2);
                val = checkFace(e, f2);
                if val
                    del.length = del.length + 1;
                    del.i(del.length) = e;
                    del.j(del.length) = f2;
                    del.val(del.length) = val;
                end
            end
        end
    end
end


function searchInVertexY2TriangleX(e)
    offset = 3 * numFacesXxNumFacesY + 6* numFacesY*numEdgesX+ numVerticesX * numFacesY;
    f = offset + (eToEdgesYTranslator(e)-1) * numFacesX;
    f1 = f + locEdgesXinFX(eToEdgesXTranslator(e), 1);
    % first orientation of the edge
    val = checkFace(e, f1);
    if val
        del.length = del.length + 1;
        del.i(del.length) = e;
        del.j(del.length) = f1;
        del.val(del.length) = val;
    end
    % second orientation of the edge
    f2 = f + locEdgesXinFX(eToEdgesXTranslator(e), 2);
    val = checkFace(e, f2);
    if val
        del.length = del.length + 1;
        del.i(del.length) = e;
        del.j(del.length) = f2;
        del.val(del.length) = val;
    end
    
end


function searchInEdgesYToTriangleX(e)
    offset = 3 * numFacesXxNumFacesY + 6* numFacesY*numEdgesX+ numVerticesX * numFacesY + ...
    numVerticesY * numFacesX;
    
    % the i loop accounts for the possible rotations
    for i = 0:3-1 
        rotPos = i * 2 * numEdgesY * numFacesX;
        
        % the ii loop accounts for different possibilities the triangles
        % can be constructet from one edge
        for ii = 0:2-1 
            edge2TriPos = ii * numEdgesY * numFacesX;
            
            f = offset + edge2TriPos + rotPos + (eToEdgesYTranslator(e)-1) * numFacesX;
            f1 = f + locEdgesXinFX(eToEdgesXTranslator(e), 1);
            val = checkFace(e, f1);
            if val
                del.length = del.length + 1;
                del.i(del.length) = e;
                del.j(del.length) = f1;
                del.val(del.length) = val;
            end
            f2 = f + locEdgesXinFX(eToEdgesXTranslator(e), 2);
            val = checkFace(e, f2);
            if val
                del.length = del.length + 1;
                del.i(del.length) = e;
                del.j(del.length) = f2;
                del.val(del.length) = val;
            end
        end
    end
end


function searchDegEdgesYInEdgesY2TriangleX(e)
    offset = 3 * numFacesXxNumFacesY + 6* numFacesY*numEdgesX+ numVerticesX * numFacesY + ...
            numVerticesY * numFacesX;
    
    % we don't know how often a vertex appears in the corresponding edge
    % matrix of a shape => this functions searches for the vertices
    [idxes, numIdxes] = findVerticesInEdgesMatrix(eToEdgesYTranslator(e), 0, 1);
    
    % the i loop accounts for the possible rotations
    for i = 0:3-1 
        rotPos = i * 2 * numEdgesY * numFacesX;
        
        % the vertex index can either be in the first or the second element
        % of an edge => k == 1 first element; k == 2 second element
        for k = 1:2
            for idx = 1:numIdxes(k)
                edge2TriPos = (idxes(idx, k)-1) * numFacesX + (k-1) ...
                    * numEdgesY * numFacesX;

                f = offset + edge2TriPos + rotPos;
                % first orientation of the edge
                f1 = f + locEdgesXinFX(eToEdgesXTranslator(e), 1);
                val = checkFace(e, f1);
                if val
                    del.length = del.length + 1;
                    del.i(del.length) = e;
                    del.j(del.length) = f1;
                    del.val(del.length) = val;
                end
                % first orientation of the edge
                f2 = f + locEdgesXinFX(eToEdgesXTranslator(e), 2);
                val = checkFace(e, f2);
                if val
                    del.length = del.length + 1;
                    del.i(del.length) = e;
                    del.j(del.length) = f2;
                    del.val(del.length) = val;
                end
            end
        end
    end
end


%% Helper functions
function val = checkFace(e, f)
    if f <= -1
        val = 0;
        return 
    end
    productFace = [FXCombo(f, :), FYCombo(f, :)] ;

    edge1 = productFace([1 2 4 5]);
    val = checkEdge(edge1, e);
    if val
        return;
    end
    edge2 = productFace([2 3 5 6]);
    val = checkEdge(edge2, e);
    if val
        return;
    end
    edge3 = productFace([3 1 6 4]);
    val = checkEdge(edge3, e);
end

function val = checkEdge(edge, e) 
    val = 0;
    if old
        if isequal(E(e, :), edge) 
            val = 1;
        elseif isequal(E(e, :), edge([2 1 4 3])) 
            val = -1;
        end
    else
        if isequal(E(e, :), edge)
            val = 1;
        elseif isequal(E(e, :), edge([1 2 4 3]))
            val = -1;
        elseif isequal(E(e, :), edge([2 1 4 3]))
            val = -1;
        elseif isequal(E(e, :), edge([2 1 3 4]))
            val = 1;
        end
    end   
end


function [idxes, numIdxes] = findVerticesInEdgesMatrix(vertexIdx, reset, shape)
%findVerticesInEdgesMatrix finds a vertex represented by its idx within the
%                           edge matrix of the specified shape
%
% This function stores its result and just if a reset call occurs or the
% vertexIdx changes, it updates its values.
%
% Inputs:
%   vertexIdx:      integer of the vertex for which we should search
%   reset:          bool flag which indicates if we should reset the
%                   function
%   shape:          bool to distiguish in which edge matrix we should
%   search
% shape == 0 => shapeX
% shape == 1 => shapeY
    if reset
        oldVertexIdx = -1;
        perIdxes = zeros(50, 2);
        perNumIdxes = zeros(1, 2);
        return;
    end
    
    % if there is no change return the old value
    if vertexIdx == oldVertexIdx
        idxes = perIdxes;
        numIdxes = perNumIdxes;
        return;
    end
    % we don't know how many edges will contain the vertex so we go with a
    % pretty high value of 50 here to avoid reallocation of array space
    idxes = zeros(50, 2); 
    numIdxes = zeros(1, 2);
    if ~shape
        EDG = edgesX; 
        nEDG = numEdgesX;
    else
        EDG = edgesY; 
        nEDG = numEdgesY;
    end
    
    % iterate through the edges of shape
    for ee = 1:nEDG
        if EDG(ee, 1) == vertexIdx
            numIdxes(1) = numIdxes(1) + 1;
            idxes(numIdxes(1), 1) = ee;
        end
        if EDG(ee, 2) == vertexIdx
            numIdxes(2) = numIdxes(2) + 1;
            idxes(numIdxes(2), 2) = ee;
        end
        
    end
    perIdxes = idxes;
    perNumIdxes = numIdxes;
end

end




