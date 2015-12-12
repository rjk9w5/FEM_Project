function [globalStiffnessMatrix, globalForceVector] = AssembleGlobalMatrices(meshNodes,meshElems,degreesOfFreedomPerNode)

[numberOfElements, nodesPerElement] = size(meshElems);
[numberOfNodes, spatialDimension] = size(meshNodes);
degreesOfFreedomPerElement = degreesOfFreedomPerNode*nodesPerElement;

globalDegreesOfFreedom = degreesOfFreedomPerNode*numberOfNodes;
globalStiffnessMatrix = zeros(globalDegreesOfFreedom, ...
                              globalDegreesOfFreedom);
globalForceVector = zeros(globalDegreesOfFreedom,1);

h = waitbar(0,'Assembling...');

for elem=1:numberOfElements;
    waitbar(elem/numberOfElements,h);
    nodes = meshElems(elem,1:nodesPerElement);
    [localStiffnessMatrix, localForceVector, B, tangentMatrix] = ...
        TriangleElementMatrices(nodes,meshNodes);
    
    for i=1:nodesPerElement
        local_i = (i-1)*degreesOfFreedomPerNode + 1;
        global_i = (nodes(i)-1)*degreesOfFreedomPerNode + 1;
        for j=1:nodesPerElement
            local_j = (j-1)*degreesOfFreedomPerNode + 1;
            global_j = (nodes(j)-1)*degreesOfFreedomPerNode + 1;
 
    globalStiffnessMatrix(global_i:(global_i+degreesOfFreedomPerNode-1), ...
                          global_j:(global_j+degreesOfFreedomPerNode-1)) ...
  = globalStiffnessMatrix(global_i:(global_i+degreesOfFreedomPerNode-1), ...
                          global_j:(global_j+degreesOfFreedomPerNode-1)) ...
    + localStiffnessMatrix(local_i:(local_i+degreesOfFreedomPerNode-1), ...
                           local_j:(local_j+degreesOfFreedomPerNode-1));
        end % loop over j
              
      globalForceVector(global_i:(global_i+degreesOfFreedomPerNode-1),1) ...
    = globalForceVector(global_i:(global_i+degreesOfFreedomPerNode-1),1) ...
      + localForceVector(local_i:(local_i+degreesOfFreedomPerNode-1),1);
    end % loop over i
    
end % loop over elements
delete(h);
end % function