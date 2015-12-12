close all

% specimen geometry
w = 4; % m
h = 1; % m
r = 1; % m
A = 2*w - 4*r;
L = 0.2*A; % m
%%declaring a tollerence

% boundary condition
appliedDisplacement= 0.001*L; % apply 0.1% cross-head strain

% simulation parameters
meshSize = 0.7; % m
deformationScaleFactorForPlotting = 100; % Multiplies displacement to 
% that small displacements are visible in plots.

% generate dog bone mesh using Dist Mesh program
[meshNodes,meshElems]=GenerateDogBoneMesh(w,h,L+2*r,r,meshSize);

% For plane stress, we just have 1 degree of freedom per node
degreesOfFreedomPerNode = 2;

% Assemble global stiffnes matrix and global force vector
[globalStiffnessMatrix, globalForceVector] = ...
    AssembleGlobalMatrices(meshNodes,meshElems,degreesOfFreedomPerNode);

% Store stiffness matrix to compute nodal forces later
initialGlobalStiffnessMatrix = globalStiffnessMatrix; 

% Apply boundary conditions to global stiffness matrix and initialize
% some parameters.
matSize = size(globalForceVector);
globalDegreesOfFreedom = matSize(1,1);
tol = 10^-4;
zeroVector = zeros(globalDegreesOfFreedom,1);
numberOfNodes = size(meshNodes);

% modify global matrices
for nodeIndex=1:numberOfNodes
    uIndex = (nodeIndex-1)*degreesOfFreedomPerNode + 1;
    vIndex = (nodeIndex-1)*degreesOfFreedomPerNode + 2;
    x = meshNodes(nodeIndex,1);
    y = meshNodes(nodeIndex,2);
    
    % TODO (first): apply displacement boundary conditions // done! :)
    % fix the displacement to be zero on the bottom
    
    if abs(y-0) <= tol %%checks to see if the node is at the bottom
        
        globalStiffnessMatrix(uIndex,:)=0;
        globalStiffnessMatrix(uIndex,uIndex)=1;
        globalForceVector(uIndex,1)=0;
        globalStiffnessMatrix(vIndex,:)=0;
        globalStiffnessMatrix(vIndex,vIndex)=1;
        globalForceVector(vIndex,1)=0;
        
    end
    
    if abs(y-(2*h+2*r+L)) <= tol %%checks to see if the node is on the top of the mesh
        
        globalStiffnessMatrix(uIndex,:)=0;
        globalStiffnessMatrix(uIndex,uIndex)=1;
        globalForceVector(uIndex,1)=0;
        globalStiffnessMatrix(vIndex,:)=0;
        globalStiffnessMatrix(vIndex,vIndex)=1;
        globalForceVector(vIndex,1)=appliedDisplacement;
        
    end
    

end

% solve for nodal displacements
U = globalStiffnessMatrix\globalForceVector;

% compute total forces on bottom
nodalForces = initialGlobalStiffnessMatrix*U; % a vector containing nodal
% forces
totalForce = 0;
for nodeIndex=1:numberOfNodes
    uIndex = (nodeIndex-1)*degreesOfFreedomPerNode + 1;
    vIndex = (nodeIndex-1)*degreesOfFreedomPerNode + 2;
    x = meshNodes(nodeIndex,1);
    y = meshNodes(nodeIndex,2);
    
    % TODO (second): sum up total nodal forces on the bottom surface Done :)
    % if the current node is on the bottom surface, add it's force in the
    % vertical direction to totalForce
    
    if abs(y)<=tol  %%checks to see if the node is on the bottom
        
        totalForce=totalForce+nodalForces(nodeIndex);
        
    end

end

% TODO (third): compute engineering stress and strain  Done :)
% compute engineerings stress and strain
engineeringStress=-totalForce/A;
engineeringStrain=appliedDisplacement/L;

% TODO (fourth): compute apparent Youngs modulus Done :)
E=engineeringStress/engineeringStrain

% compute displacement and strain field for plotting
[stresses,strains] = NodalStressStrain(meshNodes,meshElems,U,degreesOfFreedomPerNode);

% extract u and v nodal displacements for plotting
u = U(1:2:globalDegreesOfFreedom,1);
v = U(2:2:globalDegreesOfFreedom,1);

% plot solution
PlotFieldonMesh(meshNodes,meshElems,u);
axis([-w w 0 2*h+L+2*r]);
title('u displacement');
xlabel('x (m)');
ylabel('y (m)');

PlotFieldonMesh(meshNodes,meshElems,v);
axis([-w w 0 2*h+L+2*r]);
title('v displacement');
xlabel('x (m)');
ylabel('y (m)');

depl = [u, v];

PlotFieldonDefoMesh(meshNodes,meshElems,deformationScaleFactorForPlotting,depl,v);
title('v displacement');
xlabel('x (m)');
ylabel('y (m)');

PlotFieldonMesh(meshNodes,meshElems,stresses(1,:));
axis([-w w 0 2*h+L+2*r]);
title('\sigma_{xx}');
xlabel('x (m)');
ylabel('y (m)');

PlotFieldonMesh(meshNodes,meshElems,stresses(2,:));
axis([-w w 0 2*h+L+2*r]);
title('\sigma_{yy}');
xlabel('x (m)');
ylabel('y (m)');

PlotFieldonMesh(meshNodes,meshElems,stresses(3,:));
axis([-w w 0 2*h+L+2*r]);
title('\sigma_{xy}');
xlabel('x (m)');
ylabel('y (m)');