addpath(genpath('distmesh'));
addpath(genpath('postprocessing'));

% generate mesh of square for patch test
meshNodes = zeros(9,2);
meshNodes(1,:) = [0,0];
meshNodes(2,:) = [0.5,0];
meshNodes(3,:) = [1,0];
meshNodes(4,:) = [0,0.5];
meshNodes(5,:) = [0.45,0.4];
meshNodes(6,:) = [1,0.5];
meshNodes(7,:) = [0,1];
meshNodes(8,:) = [0.5,1];
meshNodes(9,:) = [1,1];

meshElems = zeros(8,3);
meshElems(1,:) = [1, 2, 4];
meshElems(2,:) = [2, 5, 4];
meshElems(3,:) = [2, 3, 5];
meshElems(4,:) = [3, 6, 5];
meshElems(5,:) = [4, 5, 7];
meshElems(6,:) = [5, 8, 7];
meshElems(7,:) = [5, 6, 8];
meshElems(8,:) = [6, 9, 8];

% applied Displacement on top
appliedDisplacement= 0.001;

% For plane stress, we just have 1 degree of freedom per node
degreesOfFreedomPerNode = 2;

% Assemble global stiffnes matrix and global force vector
[globalStiffnessMatrix, globalForceVector] = ...
    AssembleGlobalMatrices(meshNodes,meshElems,degreesOfFreedomPerNode);

% Apply boundary conditions
% initialize some parameters
matSize = size(globalForceVector);
globalDegreesOfFreedom = matSize(1,1);
tol = 10^-6;
zeroVector = zeros(globalDegreesOfFreedom,1);
numberOfNodes = size(meshNodes);

% modify global matrices
for nodeIndex=1:numberOfNodes
    uIndex = (nodeIndex-1)*degreesOfFreedomPerNode + 1;
    vIndex = (nodeIndex-1)*degreesOfFreedomPerNode + 2;
    x = meshNodes(nodeIndex,1);
    y = meshNodes(nodeIndex,2);
    
    % fix the displacement to be zero on the bottom
    if abs(y) < tol
        if abs(x) < tol
        globalStiffnessMatrix(uIndex,:) = zeroVector;
        globalStiffnessMatrix(uIndex,uIndex) = 1;
        globalForceVector(uIndex,1) = 0;
        end
        globalStiffnessMatrix(vIndex,:) = zeroVector;
        globalStiffnessMatrix(vIndex,vIndex) = 1;
        globalForceVector(vIndex,1) = 0;
    end
    % apply the displacement at the top
    if abs(y-1) < tol
        globalStiffnessMatrix(vIndex,:) = zeroVector;
        globalStiffnessMatrix(vIndex,vIndex) = 1;
        globalForceVector(vIndex,1) = appliedDisplacement;
    end

end

% solve for nodal displacements
U = globalStiffnessMatrix\globalForceVector;

% compute displacement and strain field for plotting
[stresses,strains] = NodalStressStrain(meshNodes,meshElems,U,degreesOfFreedomPerNode);

% extract u and v nodal displacements for plotting
u = U(1:2:globalDegreesOfFreedom,1);
v = U(2:2:globalDegreesOfFreedom,1);

% plot solution
PlotFieldonMesh(meshNodes,meshElems,u);
title('u displacement');
xlabel('x (m)');
ylabel('y (m)');

PlotFieldonMesh(meshNodes,meshElems,v);
title('v displacement');
xlabel('x (m)');
ylabel('y (m)');

strains