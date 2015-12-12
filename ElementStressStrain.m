function [elementStress, elementStrain] = ElementStressStrain(nodes,meshNodes,elementDisplacements)
  
[localStiffnessMatrix, localForceVector, B, tangentMatrix] = ...
        TriangleElementMatrices(nodes,meshNodes);

% TODO: compute element stress and strain (replace next two lines)
elementStrain = ones(3,1);
elementStress = ones(3,1);

end