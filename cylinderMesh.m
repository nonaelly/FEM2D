function [nodes, elements4Node, nodeBou, elemBou] = cylinderMesh(rInner, rOuter, numR, numTheta)
% Generate a 2D mesh for a 1/4 thick-walled cylinder and extract boundary nodes
% rInner: Inner radius
% rOuter: Outer radius
% numR: Number of nodes in the radial direction
% numTheta: Number of nodes in the angular direction

% Create radial and angular grid points
r = linspace(rInner, rOuter, numR);
theta = linspace(0, pi/2, numTheta);

% Create structured mesh grid (polar coordinates to Cartesian)
[R, Theta] = meshgrid(r, theta);
X = R .* cos(Theta);
Y = R .* sin(Theta);
Z = zeros(size(X));

% Plot the mesh
figure;
plot(X, Y, 'k'); % Radial lines
hold on;
plot(X', Y', 'k'); % Angular lines
axis equal;
title('1/4 Thick-Walled Cylinder Structured Mesh');

% Generate node numbering
nodes = [(1:numel(X))', X(:), Y(:), Z(:)]; % Convert the grid into a list of nodes

% Generate 4-node quadrilateral elements
k = 1;
elements4Node = [];
elemBou = cell(4, 1);
for i = 1:numTheta-1
    for j = 1:numR-1
        n1 = i + (j-1) * numTheta;         % Bottom-left
        n2 = i + j * numTheta;     % Bottom-right
        n3 = (i+1) + j * numTheta;         % Top-right
        n4 = (i+1) + (j-1) * numTheta;             % Top-left
        % 检查单元是否在边界上，并记录对应的 elemBou
        if j == 1 % Inner boundary
            elemBou{1} = [elemBou{1}; k, 1]; % 单元编号，边界面（内侧）
        elseif j == numR-1 % Outer boundary
            elemBou{3} = [elemBou{3}; k, 3]; % 单元编号，边界面（外侧）
        end

        % theta = 0 边界
        if i == 1
            elemBou{2} = [elemBou{2}; k, 2]; % 单元编号，边界面（theta=0）
        end

        % theta = pi/2 边界
        if i == numTheta-1
            elemBou{4} = [elemBou{4}; k, 4]; % 单元编号，边界面（theta=pi/2）
        end

        elements4Node = [elements4Node; k, 1, n1, n2, n3, n4];
        k = k + 1;
    end
end

% Extract boundary nodes
innerBoundary = nodes(1:numTheta, 1:3);       % Inner radius boundary
outerBoundary = nodes((numTheta*(numR-1)+1):end, 1:3); % Outer radius boundary
theta0Boundary = nodes(1:numTheta:end, 1:3); % theta = 0 boundary
theta90Boundary = nodes(numTheta:numTheta:end, 1:3);   % theta = pi/2 boundary

nodeBou = cell(4, 1);
nodeBou{1} = innerBoundary;
nodeBou{2} = theta0Boundary;
nodeBou{3} = outerBoundary;
nodeBou{4} = theta90Boundary;


% Display boundaries
%     disp('Inner Boundary Nodes:');
%     disp(innerBoundary);
%     disp('Outer Boundary Nodes:');
%     disp(outerBoundary);
%     disp('Theta = 0 Boundary Nodes:');
%     disp(theta0Boundary);
%     disp('Theta = pi/2 Boundary Nodes:');
%     disp(theta90Boundary);
%
%     % Display element connectivity
%     disp('4-Node Quadrilateral Elements:');
%     disp(elements4Node);
end
