% Code to reproduce Fig 4.4(a) using Matlab PDE Toolbox - Mesh Generating
% function to create Mesh in Fig 4.3(c)

clear all
close all

bndryN = 16;
intrN = 48;
totalN = bndryN + intrN;

% % =======================================================================
% % =======================================================================
% Determine the Test Points
r = 0.3;
[p,e,t] = initmesh('circleg','Hmax',r); % create Circular Mesh with 64 points
tstX = p(1,:)';
tstY = p(2,:)';
% % =======================================================================
% Determine the Points of Grid
r = 0.35;
[p,e,t] = initmesh('circleg','Hmax',r); % create Circular Mesh with 64 points

% p : DESIRED CIRCULAR MESH for RBF INTERPOLATION.
gridX = p(1,:)';
gridY = p(2,:)';

% ------------------------------------------------------
% plot Unit Circle
theta = 0:0.01:2*pi;
x = cos(theta);
y = sin(theta);
fig1 = figure(1)
plot(x,y)
hold on
% plot the mesh-points
scatter(gridX, gridY,'r')
hold off

% % =======================================================================
    % % =======================================================================
    %[g f] = evalPDEatPoints(p);
    % --------------------------------------------
    % g = u(x,y) =  A ./ (A + (x - a).^2 + b* y.^2);
    % f = Laplacian[u(x,y)] = -2*A*(A + (x - a).^2 + b* y.^2).^(-3)
    % .*(A(b+1) + (x-a)^2(b-3) + y^2(b -3b^2))
    % --------------------------------------------
    % Boundary Constraint
    g = @(X, Y) (100 ./(100 + (X - 0.2).^2 + 2* Y.^2) );
    % Interior Constraint
    f = @(X,Y) (-2*100*(100*(2+1) + (2-3)*(X-0.2).^2 + (2 - 3*2^2)*Y.^2).*(100 + (X - 0.2).^2 + 2*Y.^2).^(-3));
    
    
% % =======================================================================
% % =======================================================================
% Solving PDE using RBF
        

epsilons = 0.001:0.001:10;
maxError = zeros(length(epsilons),3);

for i = 1:3
    switch i
        case 1
            rbf = 'GA'
        case 2
            rbf = 'IQ'
        case 3 
            rbf = 'MQ'
    end

for count = 1:length(epsilons)
    ep = epsilons(count);
   
    % % =======================================================================
    % % =======================================================================
    % Kansa's Formulation
    % --------------------------------------------------
    % Determine Vectors
    
    gridR = gridX.^2 + gridY.^2;
    
    % Identify index of points that are located on UNIT CIRCLE BOUNDARY
    bdryID = find(gridR==1);
    intrID = find(gridR<1);
    
    % Rearrange Points such that first values are on Boundary and remaining
    % are in the interior
    
    newOrder = [bdryID; intrID];
    gridX = gridX(newOrder);
    gridY = gridY(newOrder); 
    
    numBoundaryPts = length(bdryID);
    
    % Evaluate Boudary Constraint
    gridX_BND = gridX(1:numBoundaryPts);
    gridY_BND = gridY(1:numBoundaryPts);
    gGrid = g(gridX_BND,gridY_BND);
    
    % Evaluate Interior Constraint
    gridX_INT = gridX(numBoundaryPts +1:end);
    gridY_INT = gridY(numBoundaryPts +1:end);
    fGrid = f(gridX_INT,gridY_INT);
        
    % ---------------------------------------------------
    % Compute R for the A matrix
    
    [x1, x2] = meshgrid(gridX);
    [y1, y2] = meshgrid(gridY);
    d2 = (x1 - x2).^2 + (y1 - y2).^2;
    r = sqrt(d2);
    
    % Generate the A matrix for points on Boundary
    A = fi(rbf,ep, r(1:numBoundaryPts,:) );
    % Generate Laplacian_A matrix for points in Interior
    LA = Lfi(rbf,ep, r(numBoundaryPts+1:end,:) );
    
    % Matrix for Kansa's Method
    Ahat = [A; LA];
    
    % Evaluate the weighting coefficients
    lambda = Ahat\[gGrid;fGrid];
    
    % Grid points tst_i - grid_j
    [tX gX] = ndgrid(tstX, gridX);
    [tY gY] = ndgrid(tstY, gridY);
    
    tstR2 = (tX-gX).^2 + (tY - gY).^2;
    tstR = sqrt(tstR2);
    
    
    tstA = fi(rbf,ep,tstR);
    
%     fprintf('tstU_rbf\n')
    tstU_rbf = tstA*lambda;
    size(tstU_rbf);
    
%     fprintf('tstU_true\n')
    tstU_true = g(tstX,tstY);
    size(tstU_true);
    
    error = tstU_true - tstU_rbf;
    maxError(count,i) = max(abs(error));
end

end
fig2 = figure(2)
plot(log10(epsilons),log10(maxError(:,1)),'--',log10(epsilons),log10(maxError(:,2)),'-.', log10(epsilons),log10(maxError(:,3)),'-' )
title('Max norm Error')
xlabel('log10(epsilon)')
ylabel('log10(error)')
legend('GA','MQ','IQ')


fig3 = figure(3)
plot(log10(epsilons),log10(maxError(:,1)),'--g')
hold on
plot(log10(epsilons),log10(maxError(:,2)),'-.b')
hold on
plot(log10(epsilons),log10(maxError(:,3)),'-r' )
hold off
title('Max norm Error')
xlabel('log10(epsilon)')
ylabel('log10(error)')
