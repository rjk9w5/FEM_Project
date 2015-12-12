function [nodalStresses, nodalStrains] = NodalStressStrain(meshNodes,meshElems,globalDisplacements,degreesOfFreedomPerNode)

[numberOfElements, nodesPerElement] = size(meshElems);
[numberOfNodes, spatialDimension] = size(meshNodes);
degreesOfFreedomPerElement = degreesOfFreedomPerNode*nodesPerElement;

nodalStresses = zeros(3,numberOfNodes); % Matrix for storing the column 
% vector of the Voigt notation stress for each node, which correspond to
% each column.
nodalStrains = zeros(3,numberOfNodes); % Matrix for storing the column 
% vector of the Voigt notation strain for each node, which correspond to
% each column.
m = zeros(1,numberOfNodes); % Row vector for storing the
% number of elements that share each node.

h = waitbar(0,'Computing nodal stresses and strains...');

% loop over all the elements
for elem=1:numberOfElements
    waitbar(elem/numberOfElements,h)
    nodes = meshElems(elem,:);
    
    % assemble the element displacement vector
    elementDisplacements = zeros(nodesPerElement*degreesOfFreedomPerNode,1);
    for i=1:nodesPerElement
        node = nodes(i);
        elementDisplacements(((i-1)*degreesOfFreedomPerNode + 1):(i*degreesOfFreedomPerNode),1) ...
            = globalDisplacements(((node-1)*degreesOfFreedomPerNode + 1):(node*degreesOfFreedomPerNode),1);

    end
    % Using the element displacements, compute the stress and strain of
    % the element.
    [elementStress, elementStrain] = ElementStressStrain(nodes,meshNodes,elementDisplacements);
    
    for i=1:nodesPerElement
       node = nodes(i);
       % TODO (first): add the stress and strain from the current element to the
       % stress and strain of elem. Also, increment the number of
       % elementsPerNode for each node on this element.
       
       nodalStresses(:,node) = nodalStresses(:,node) + elementStress(:);
       nodalStrains(:, node) = nodalStrains(:,node) + elementStrain(:);
       m(1,node) = m(1,node) + 1;

    end
end
delete(h);

for i=1:numberOfNodes
   % TODO (second): compute the averate stress and strain at each node by diving the
   % sum by the number of adjacent elements. Replace the next two lines.
   nodalStresses(:,i) = ones(3,1);
   nodalStrains(:,i) = ones(3,1);
end

end