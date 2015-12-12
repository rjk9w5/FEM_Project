function [elementStress, elementStrain] = ElementStressStrain(nodes,meshNodes,elementDisplacements)
  
[localStiffnessMatrix, localForceVector, B, tangentMatrix] = ...
        TriangleElementMatrices(nodes,meshNodes);

% TODO: compute element stress and strain (replace next two lines)
elementStrain = B*elementDisplacements;
elementStress = tangentMatrix*elementStrain;

end