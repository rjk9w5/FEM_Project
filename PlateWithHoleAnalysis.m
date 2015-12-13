close all
clear

addpath(genpath('distMesh'));
addpath(genpath('postProcessing'));

% input mesh parameters
sideLength = 2;
radiusOfHole = 0.2;
holeMeshSize = 0.2*radiusOfHole;
maximumMeshSize = sideLength/10;

% input boundary condition
appliedTensileStress = 10^8; % Pa
displacementFactor = 200;

% generate mesh using Dist Mesh
fd=@(p) ddiff(drectangle(p,-sideLength/2,sideLength/2,-sideLength/2,sideLength/2),dcircle(p,0,0,radiusOfHole));
fh=@(p) holeMeshSize+dcircle(p,0,0,radiusOfHole)*(maximumMeshSize-holeMeshSize)/(1-radiusOfHole);
[meshNodes,meshElems]=distmesh2d(fd,fh,holeMeshSize,[-sideLength/2,-sideLength/2;sideLength/2,sideLength/2],[-1,-1;-1,1;1,-1;1,1]);

% Compute numder of nodes on top and bottom (used for computing applied
% nodal forces).
numberOfNodes = size(meshNodes);
tol=10^-6;
numNodesOnTop = 0;
numNodesOnBottom = 0;
for nodeIndex=1:numberOfNodes
   y = meshNodes(nodeIndex,2);
   
   % TODO (first): compute number of nodes on top and bottom surfaces
   % See below for (wrong) example. Replace the following code.
   if abs(y-1) < tol
       numNodesOnTop = numNodesOnTop + 1;
   end

   if abs(y+1) < tol
       numNodesOnBottom = numNodesOnBottom + 1;
   end

end

numElemsOnTop = numNodesOnTop - 1;
numElemsOnBottom = numNodesOnBottom - 1;
forcePerElementTop = appliedTensileStress*sideLength/numElemsOnTop;
forcePerElementBottom = appliedTensileStress*sideLength/numElemsOnBottom;

% For temperature, we just have 1 degree of freedom per node
degreesOfFreedomPerNode = 2;

% Assemble global stiffnes matrix and global force vector
[globalStiffnessMatrix, globalForceVector] = ...
    AssembleGlobalMatrices(meshNodes,meshElems,degreesOfFreedomPerNode);

initialGlobalStiffnessMatrix = globalStiffnessMatrix; % Store initial 
% stiffness matrix for computing nodal forces later.

% Apply boundary conditions
% initialize some parameters
matSize = size(globalForceVector);
globalDegreesOfFreedom = matSize(1,1);
tol = 10^-6;
zeroVector = zeros(globalDegreesOfFreedom,1);
numberOfNodes = size(meshNodes);

% modify global matrices
for nodeIndex=1:numberOfNodes
    xForceIndex = (nodeIndex-1)*degreesOfFreedomPerNode + 1;
    yForceIndex = (nodeIndex-1)*degreesOfFreedomPerNode + 2;
    x = meshNodes(nodeIndex,1);
    y = meshNodes(nodeIndex,2);
    r = sqrt(x^2+y^2);
    
    % TODO (second): apply force boundary and displcement boundary conditions
    % See below for (wrong) example. Need to replace.
    if (abs(y + 1) < tol )
      % apply half since in edge
      if(abs(x + 1) < tol || abs(x - 1) < tol)
        globalForceVector(yForceIndex,1) = globalForceVector(yForceIndex,1) - forcePerElementTop/2;
      else
        globalForceVector(yForceIndex,1) = globalForceVector(yForceIndex,1) - forcePerElementTop;
      end
    elseif(abs(y - 1) < tol)
      if(abs(x + 1) < tol || abs(x - 1) < tol)
        globalForceVector(yForceIndex,1) = globalForceVector(yForceIndex,1) + forcePerElementTop/2;
      else
        globalForceVector(yForceIndex,1) = globalForceVector(yForceIndex,1) + forcePerElementTop;
      end
    end

    if(abs(x-1) < tol && abs(y+1) < tol)
      globalStiffnessMatrix(yForceIndex,:) = 0;
      globalStiffnessMatrix(yForceIndex,yForceIndex) = 1;
      globalForceVector(yForceIndex,1) = 0;
    end
    if(abs(x+1) < tol && abs(y+1) < tol)
      globalStiffnessMatrix(yForceIndex,:) = 0;
      globalStiffnessMatrix(yForceIndex,yForceIndex) = 1;
      globalForceVector(yForceIndex,1) = 0;

      globalStiffnessMatrix(xForceIndex,:) = 0;
      globalStiffnessMatrix(xForceIndex,xForceIndex) = 1;
      globalForceVector(xForceIndex,1) = 0;
    end

end

% solve for nodal temperatures
U = globalStiffnessMatrix\globalForceVector;

% compute stress and strain from solution
[stresses,strains] = NodalStressStrain(meshNodes,meshElems,U,degreesOfFreedomPerNode);

% computer error of stress intensity factor
relativeHoleMeshSize = holeMeshSize/radiusOfHole
relativePercentErrorOfStress = 100*abs(max(stresses(2,:)) - 3*appliedTensileStress)/(3*appliedTensileStress)

u = U(1:2:globalDegreesOfFreedom,1);
v = U(2:2:globalDegreesOfFreedom,1);

% plot solution
PlotFieldonMesh(meshNodes,meshElems,u);
axis([-1 1 -1 1]);
title('u displacement');
xlabel('x (m)');
ylabel('y (m)');

PlotFieldonMesh(meshNodes,meshElems,v);
axis([-1 1 -1 1]);
title('v displacement');
xlabel('x (m)');
ylabel('y (m)');

depl = [u, v];

PlotFieldonDefoMesh(meshNodes,meshElems,displacementFactor,depl,v);
title('v displacement');
xlabel('x (m)');
ylabel('y (m)');

PlotFieldonMesh(meshNodes,meshElems,stresses(1,:));
axis([-1 1 -1 1]);
title('\sigma_{xx}');
xlabel('x (m)');
ylabel('y (m)');

PlotFieldonMesh(meshNodes,meshElems,stresses(2,:));
axis([-1 1 -1 1]);
title('\sigma_{yy}');
xlabel('x (m)');
ylabel('y (m)');

PlotFieldonMesh(meshNodes,meshElems,stresses(3,:));
axis([-1 1 -1 1]);
title('\sigma_{xy}');
xlabel('x (m)');
ylabel('y (m)');