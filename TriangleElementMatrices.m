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
x(1) = nodeMesh(nodes(1),1);
y(1) = nodeMesh(nodes(1),2);
x(2) = nodeMesh(nodes(2),1);
y(2) = nodeMesh(nodes(2),2);
x(3) = nodeMesh(nodes(3),1);
y(3) = nodeMesh(nodes(3),2);

derivativeMat = [-1 1 0; -1 0 1];

% TODO (first): compute the tangent matrix

tangentMatrix = YoungModulus / (1 - PoissonRatio) * [1 PoissonRatio 0; PoissonRatio 1 0; 0 0 (1 - PoissonRatio)/2];
% Replace next line
%tangentMatrix = ones(3,3);

% compute the local stiffness matrix
localStiffnessMatrix = zeros(numberOfNodesPerElement*degreesOfFreedomPerNode, ...
                             numberOfNodesPerElement*degreesOfFreedomPerNode);
for qp = 1:numberOfQuadraturePoints
    xi_qp = quadraturePoints(qp,1);
    eta_qp = quadraturePoints(qp,2);
    
    % TODO (second): compute the B matrix and Jacobian, J (replace next two
    % lines.

    J = [dot(x,derivativeMat(1,:)) dot(y,derivativeMat(:,1)); 
        dot(x,derivativeMat(2,:)), dot(y,derivativeMat(2,:))];
    derivativeMatxy = J\derivativeMat;
    B = zeros(3,6);
    B(1, 1:2:end) = derivativeMatxy(1,:);
    B(2, 2:2:end) = derivativeMatxy(2,:);
    B(3, 1:2:end) = derivativeMatxy(2,:);
    B(3, 2:2:end) = derivativeMatxy(1,:);
    %J = ones(2,2);
    
    % Keep the line below, this is right.
    localStiffnessMatrix = localStiffnessMatrix ...
     + ((B') * tangentMatrix * B) .* det(J) .* quadratureWeigts(qp,1);
end

% assume no body force so the local force vector is zero
localForceVector = [0; 0; 0; 0; 0; 0];
                       
end