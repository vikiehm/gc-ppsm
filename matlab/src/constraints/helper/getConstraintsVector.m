function constraintsVector = getConstraintsVector(Va, Fa, Vb, Fb, degenerate, closed, old)
%getConstraintsVector constructs the RHS of the constraints

numFacesA = size(Fa, 1);
numFacesB = size(Fb, 1);
if closed
    numEdgesA = 3/2 * numFacesA; 
    numEdgesB = 3/2 * numFacesB;
else
    numEdgesA = size(getEdges(Fa, closed),1);
    numEdgesB = size(getEdges(Fa, closed),1);

end
numVerticesA = size(Va, 1);
numVerticesB = size(Vb, 1);

numProductEdges = numEdgesA * numEdgesB;   
if degenerate
    numProductEdges = numProductEdges + numEdgesA * numVerticesB ...
                    + numEdgesB * numVerticesA;
end
if old
    numProductEdges = numProductEdges + numEdgesA * numEdgesB;
end

constraintsVector = [zeros(numProductEdges, 1); ones(numFacesA + numFacesB, 1)];
end