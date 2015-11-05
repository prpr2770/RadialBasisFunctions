% Code for calculating the differentiation matrix based on RBF-FD for 
% the 'solid body rotation' convection problem described by equation (2.8) 
% in the text book.

% The governing PDE is there expressed in spherical coordinates, which
% become singular at the poles. These singularities vanish automatically
% once the problem is reformulated into (x,y,z)-based RBF form.

% Assumptions:
% 1. A unit sized sphere (a = 1) 
% 2. Rotation speed that corresponds to u_0 = 1, i.e. to a rotation time of  2 pi.

clear all; close all;
ep = 2;         % Specify shape parameter to use
al = 0.4;       % Specify rotation axis relative to spherical coordinate
                % system's axis
      
fi     = @(r) 1./(1+(ep*r).^2);             % Define the RBF (here IQ)
dfi_dr = @(r) -2*r*ep^2./(1+(ep*r).^2).^2;  % Radial derivative of the RBF

xyz = load('me_0900.txt');      % Read in the node set for the sphere; here
                                % an ME set with 900 nodes
xyz(:,1:3) = fliplr(xyz(:,1:3));% Flip to avoid a node at the north pole
N = length(xyz(:,1));           % Get the total number of nodes
warning('off','MATLAB:divideByZero')    % Code becomes simpler if we do not 
                                        % need to worry about zero divides
% Express all the node points in spherical coordinates
xyz(:,9) = xyz(:,3);                    % sin(th)
xyz(:,8) = sqrt(1-xyz(:,3).^2);         % cos(th)
xyz(:,6) = xyz(:,1)./xyz(:,8);          % cos(fi)
xyz(:,7) = xyz(:,2)./xyz(:,8);          % sin(fi)

xyz(:,4) = atan2(xyz(:,2),xyz(:,1));    % fi
xyz(:,5) = asin (xyz(:,3));             % th

% % % % Compute the Distance Table for all pair-wise distances
% % % [x1,x2] = meshgrid(xyz(:,1));           % Calculate a distance table for
% % % [y1,y2] = meshgrid(xyz(:,2));           % all pairwise distances between
% % % [z1,z2] = meshgrid(xyz(:,3));           % nodes (as measured straight
% % % d2 = (x1-x2).^2+(y1-y2).^2+(z1-z2).^2;  % through the sphere
% % % 
% % % % Determine the standard RBF-A Matrix
% % % A  = fi(sqrt(d2));  % Form the standard RBF matrix A
% % % la = inv(A);        % Obtain the exp. coeff. for the cardinal data cases

% RBF-FD Differentiation Matrix
global D;       
D = zeros(N,N);
D = sparse(D);

% Determine n-nearest neighbors
n = 20;
IDX = knnsearch(xyz(:,1:3),xyz(:,1:3),'K',n,'Distance','euclidean');

for i = 1:N % Evaluate L at xyz(i,:) when RBF centered at xyz(j,:), j=1:m
            % Note that we only need a single loop
            
    nbrs = IDX(i,:);
    
    % Compute values only to particular neighbors
    dx_dfi = -xyz(i,7)*xyz(i,8);        % To obtain the weights (which form 
    dy_dfi =  xyz(i,6)*xyz(i,8);        % the rows of the DM) we need  
    dz_dfi =  0;                        % for each RBF to calculate its
                                        % derivatives at node i with
    dx_dth = -xyz(i,6)*xyz(i,9);        % respect to fi and th.
    dy_dth = -xyz(i,7)*xyz(i,9);        % We obtain these via the chain
    dz_dth =  xyz(i,8);                 % rule after first calculating 
                                        % derivatives of the mapping from
    
                                        
    % Analysis on only the neighbors : IS THIS NEEDED HERE or globally?
    XYZ = xyz(nbrs,1:3);    % Note: XYZ(1,:) = xyz(i,1:3)
    ONE = ones(length(nbrs),1);
    size(XYZ)
    size(ONE)
    % Compute the Distance Table for all pair-wise distances
    [x1,x2] = meshgrid(XYZ(:,1));           % Calculate a distance table for
    [y1,y2] = meshgrid(XYZ(:,2));           % all pairwise distances between
    [z1,z2] = meshgrid(XYZ(:,3));           % nodes (as measured straight
    d2 = (x1-x2).^2+(y1-y2).^2+(z1-z2).^2;  % through the sphere
    R = sqrt(d2);
    Ahat = fi(R);
    
    % Computing L_h: for the Local Stencil
    r = R(1,:)';
    size(r)
    
    dh_dr = dfi_dr(r);
    size(dh_dr)
    
    dh_dx = dh_dr.*(XYZ(:,1)-XYZ(1,1)*ONE)./r;
    dh_dy = dh_dr.*(XYZ(:,2)-XYZ(1,2)*ONE)./r;
    dh_dz = dh_dr.*(XYZ(:,3)-XYZ(1,3)*ONE)./r;
    dh_dx(1) = 0; dh_dy(1) = 0; dh_dz(1) = 0;   % Reset error from divByZero

    dh_dfi = dh_dx*dx_dfi + dh_dy*dy_dfi + dh_dz*dz_dfi;
    dh_dth = dh_dx*dx_dth + dh_dy*dy_dth + dh_dz*dz_dth;
        
    L_h = -(cos(al)-tan(xyz(i,5))*xyz(i,7)*sin(al))*dh_dfi + ...
            sin(al)*xyz(i,6)*dh_dth;    % L-operator at node i evaluated 
                                        % for the different RBFs
    
    % Use the RBF-FD code here to determine the weights!   
    size(Ahat)
    size(L_h)
    stencilWts = Ahat\L_h;
    
    % Insert values into Sparse Matrix: Is this correct?
    D(i,nbrs) = stencilWts';                   % Row i of the DM computed
end

fig1 = figure(1)
E = full(D);
subplot(2,1,1)
plot(eig(E),'k.'); axis square          % Plot the eigenvalues of the DM
title('Eigenvalues: RBF-FD')
subplot(2,1,2)
spy(D)
title('Sparsity Pattern: RBF-FD')
fprintf('Check D for symmetry\n')
issymmetric(D)

fig2 = figure(2)
rcm = symrcm(D);
Drcm = D(rcm,rcm);  % Sparse matrix
Ercm = full(Drcm);  % Full matrix
subplot(2,1,1)
plot(eig(Ercm),'k.'); axis square          % Plot the eigenvalues of the DM
title('Eigenvalues: RCMK ordering ')
subplot(2,1,2)
spy(Drcm)
title('Sparsity Pattern: Reverse Cuthill-McKee Ordering ')

fprintf('Check Drcm for symmetry\n')
issymmetric(Drcm)



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% PDE - TIME STEPPING
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Convert from Cartesian Coordinates to Spherical Cordinates
ptsXYZ = xyz(:,1:3);
ptsX = xyz(:,1);
ptsY = xyz(:,2);
ptsZ = xyz(:,3);

Phi = xyz(:,4);
The = xyz(:,5);

% -------------------------------------------------------------------------
% Plot the Cosine Bell Function
m = 0:0.01:1;
R = 1/3;
h0 = 1;
ratio = m/R;
indx = (ratio<1);
hfun  = h0/2*(1+cos(pi.*ratio)).*indx;
fig3 = figure(3)
plot(m,hfun)
xlabel('r')
ylabel('Cosine Bell: h(r)')

% Evaluate the Cosine Bell for Initial Condition
r = acos(cos(The).*cos(Phi));
ratio = r/R;
indx = (ratio<1);
initH  = h0/2*(1+cos(pi.*ratio)).*indx;         % Initialize Function

% -------------------------------------------------------------------------
% Plotting using RegularizeData3D (Matlab Central)
theta = -2:0.1:2;
phi = -3.5:0.1:3.5;
Smoothness = 0.00005;

fig4 = figure(4)
z0 = RegularizeData3D(The,Phi,initH, theta, phi, 'interp', 'bicubic', 'smoothness', Smoothness);
surf(theta, phi, z0, 'facealpha', 0.50);
hold on
scatter3(The,Phi,initH, 'fill');
hold off
xlabel('\theta')
ylabel('\phi')
zlabel('h(\theta,\phi,t = 0)')
title('Initialization of Cosine Bell Curve')


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Time-Step using RK4
fractionOfRevolution = 1000;            % Ensure divisible by 4
T1rev = 2*pi;
t0 = 0;
tf = (10)*T1rev; % around 100 revolutions
dT = T1rev/fractionOfRevolution;
[t,H] = rk4_hw4(@fun_dudt_hw4, t0:dT:tf, initH, ptsXYZ);  

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Plots of Time-Revolution: RegularizeData3D (Matlab Central FileID: #46223)

fig10 = figure(10)
subplot(2,2,1)
data = initH;
z0 = RegularizeData3D(The,Phi,data, theta, phi, 'interp', 'bicubic', 'smoothness', Smoothness);
surf(theta, phi, z0, 'facealpha', 0.50);
hold on
scatter3(The,Phi,data, 'fill');
hold off
xlabel('\theta')
ylabel('\phi')
zlabel('h(\theta,\phi,t)')
title('Initialization of Cosine Bell Curve')
% -------------------------------
subplot(2,2,2)
data = H(:,fractionOfRevolution/4);
zHalfT = RegularizeData3D(The,Phi,data, theta, phi, 'interp', 'bicubic', 'smoothness', Smoothness);
surf(theta, phi, zHalfT, 'facealpha', 0.50);
hold on
scatter3(The,Phi,data,'fill')
hold off
xlabel('\theta')
ylabel('\phi')
zlabel('h(\theta,\phi,t)')
title('PDE Solution @t = T/4')
% -------------------------------
subplot(2,2,3)
data = H(:,2*fractionOfRevolution/4);
zHalfT = RegularizeData3D(The,Phi,data, theta, phi, 'interp', 'bicubic', 'smoothness', Smoothness);
surf(theta, phi, zHalfT, 'facealpha', 0.50);
hold on
scatter3(The,Phi,data,'fill')
hold off
xlabel('\theta')
ylabel('\phi')
zlabel('h(\theta,\phi,t)')
title('PDE Solution @t = T/2')
% -------------------------------
subplot(2,2,4)
data = H(:,3*fractionOfRevolution/4);
zHalfT = RegularizeData3D(The,Phi,data, theta, phi, 'interp', 'bicubic', 'smoothness', Smoothness);
surf(theta, phi, zHalfT, 'facealpha', 0.50);
hold on
scatter3(The,Phi,data,'fill')
hold off
xlabel('\theta')
ylabel('\phi')
zlabel('h(\theta,\phi,t)')
title('PDE Solution @t = 3*T/4')
% -------------------------------
% Plot Error after full revolution
fig11 = figure(11)
error = abs(H(:,fractionOfRevolution) - initH);
data = error;
Smoothness = 0.00005;

zFullT = RegularizeData3D(The,Phi,data, theta, phi, 'interp', 'bicubic', 'smoothness', Smoothness);
surf(theta, phi, zFullT, 'facealpha', 0.75);
hold on
scatter3(The,Phi,data,'fill')
hold off
xlabel('\theta')
ylabel('\phi')
zlabel('error(\theta,\phi,t)')
title('Error after 1 revolutions')

% -------------------------------
% Plot Error after 10 full revolution
fig12 = figure(12)
error = abs(H(:,end) - initH);
data = error;
Smoothness = 0.00005;

zFullT = RegularizeData3D(The,Phi,data, theta, phi, 'interp', 'bicubic', 'smoothness', Smoothness);
surf(theta, phi, zFullT, 'facealpha', 0.75);
hold on
scatter3(The,Phi,data,'fill')
hold off
xlabel('\theta')
ylabel('\phi')
zlabel('error(\theta,\phi,t)')
title('Error after 10 revolutions')


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Plotting local stencil weights
close all;
selNodes = 4;
nodes = floor(rand(selNodes,1)*N);

nbrNodes = IDX(nodes,:);
minMarkerSize = 10;
maxMarkerSize = 80;


fig30 = figure(30)

for i=1:selNodes/2
    % ------------------------------------------------
    subplot(selNodes/2,2,2*i-1)
    node = 2*i-1;
    stencilFocus = nodes(node);
    nbrs = nbrNodes(node,:);
    nbrsX = xyz(nbrs,1);
    nbrsY = xyz(nbrs,2);
    nbrsZ = xyz(nbrs,3);
    nbrWeights = D(stencilFocus,nbrs)';
    % scatter all points on sphere
    scatter(xyz(:,1),xyz(:,2),'g+')
    hold on
    % -----------------------
    pz = abs(nbrWeights);
    markerSizes = minMarkerSize + floor(pz/max(pz)*maxMarkerSize);
    red = (nbrWeights<0);
    green = zeros(size(nbrWeights));
    blue = (nbrWeights>=0);
    colorspec = [red green blue];
    % scatter Local Stencil points
    scatter(nbrsX,nbrsY,markerSizes,colorspec,'filled')
    hold off
    % ------------------------------------------------
    subplot(selNodes/2,2,2*i)
    node = 2*i;
    stencilFocus = nodes(node);
    nbrs = nbrNodes(node,:);
    nbrsX = xyz(nbrs,1);
    nbrsY = xyz(nbrs,2);
    nbrsZ = xyz(nbrs,3);
    nbrWeights = D(stencilFocus,nbrs)';
    % scatter all points on sphere
    scatter(xyz(:,1),xyz(:,2),'g+')
    hold on
    % -----------------------
    pz = abs(nbrWeights);
    markerSizes = minMarkerSize + floor(pz/max(pz)*maxMarkerSize);
    red = (nbrWeights<0);
    green = zeros(size(nbrWeights));
    blue = (nbrWeights>=0);
    colorspec = [red green blue];
    % scatter Local Stencil points
    scatter(nbrsX,nbrsY,markerSizes,colorspec,'filled')
    hold off
end
