function [localStiffnessMatrix, localForceVector, B, tangentMatrix] = TriangleElementMatrices(nodes,nodeMesh)
%%mitch was here
% element properties (e.g. aluminum)
PoissonRatio = 0.32;
YoungModulus = 70*10^9;
  
numberOfNodesPerElement = 3;
degreesOfFreedomPerNode = 2;
spatialDimension = 2;

numberOfQuadraturePoints = 1;
quadraturePoints = [1/3, 1/3];
quadratureWeigts = [1/2];

% extract nodal coordinates from mesh data
xa = nodeMesh(nodes(1),1);
ya = nodeMesh(nodes(1),2);
xb = nodeMesh(nodes(2),1);
yb = nodeMesh(nodes(2),2);
xc = nodeMesh(nodes(3),1);
yc = nodeMesh(nodes(3),2);

% TODO (first): compute the tangent matrix

tangentMatrix = E/(1-PoissonRatio) [1 PoissonRatio 0; PoissonRatio 1 0; 0 0 (1-PoissonRatio)/2];
% Replace next line
tangentMatrix = ones(3,3);

% compute the local stiffness matrix
localStiffnessMatrix = zeros(numberOfNodesPerElement*degreesOfFreedomPerNode, ...
                             numberOfNodesPerElement*degreesOfFreedomPerNode);
for qp = 1:numberOfQuadraturePoints
    xi_qp = quadraturePoints(qp,1);
    eta_qp = quadraturePoints(qp,2);
    
    % TODO (second): compute the B matrix and Jacobian, J (replace next two
    % lines.
    B = ones(3,6);
    J = ones(2,2);
    
    % Keep the line below, this is right.
    localStiffnessMatrix = localStiffnessMatrix ...
     + ((B')*tangentMatrix*B).*det(J).*quadratureWeigts(qp,1);
end

% assume no body force so the local force vector is zero
localForceVector = [0; 0; 0; 0; 0; 0];
                       
end