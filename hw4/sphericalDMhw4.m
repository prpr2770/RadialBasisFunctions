% Code for calculating the differentiation matrix based on global RBFs for 
% the 'solid body rotation' convection problem described by equation (2.8) 
% in the text book.

% The governing PDE is there expressed in spherical coordinates, which
% become singular at the poles. These singularities vanish automatically
% once the problem is reformulated into (x,y,z)-based RBF form.

% We assume a unit sized sphere (a = 1) and a rotation speed that 
% corresponds to u_0 = 1, i.e. to a rotation time of  2 pi.

clear all; close all;
ep = 2;         % Specify shape parameter to use
al = 0.4;       % Specify rotation axis relative to spherical coordinate
                % system's axis
      
fi     = @(r) 1./(1+(ep*r).^2);             % Define the RBF (here IQ)
dfi_dr = @(r) -2*r*ep^2./(1+(ep*r).^2).^2;  % Radial derivative of the RBF

xyz = load('me_0900.txt');      % Read in the node set for the sphere; here
                                % an ME set with 900 nodes
xyz(:,1:3) = fliplr(xyz(:,1:3));% Flip to avoid a node at the north pole
n = length(xyz(:,1));           % Get the total number of nodes
warning('off','MATLAB:divideByZero')    % Code becomes simpler if we do not 
                                        % need to worry about zero divides
% Express all the node points in spherical coordinates
xyz(:,9) = xyz(:,3);                    % sin(th)
xyz(:,8) = sqrt(1-xyz(:,3).^2);         % cos(th)
xyz(:,6) = xyz(:,1)./xyz(:,8);          % cos(fi)
xyz(:,7) = xyz(:,2)./xyz(:,8);          % sin(fi)

xyz(:,4) = atan2(xyz(:,2),xyz(:,1));    % fi
xyz(:,5) = asin (xyz(:,3));             % th

global D;
D = zeros(n,n);                         % Creare array to store the DM

[x1,x2] = meshgrid(xyz(:,1));           % Calculate a distance table for
[y1,y2] = meshgrid(xyz(:,2));           % all pairwise distances between
[z1,z2] = meshgrid(xyz(:,3));           % nodes (as measured straight
d2 = (x1-x2).^2+(y1-y2).^2+(z1-z2).^2;  % through the sphere
A  = fi(sqrt(d2));  % Form the standard RBF matrix A
la = inv(A);        % Obtain the exp. coeff. for the cardinal data cases

for i = 1:n % Evaluate L at xyz(i,:) when RBF centered at xyz(j,:), j=1:m
            % Note that we only need a single loop
    dx_dfi = -xyz(i,7)*xyz(i,8);        % To obtain the weights (which form 
    dy_dfi =  xyz(i,6)*xyz(i,8);        % the rows of the DM) we need  
    dz_dfi =  0;                        % for each RBF to calculate its
                                        % derivatives at node i with
    dx_dth = -xyz(i,6)*xyz(i,9);        % respect to fi and th.
    dy_dth = -xyz(i,7)*xyz(i,9);        % We obtain these via the chain
    dz_dth =  xyz(i,8);                 % rule after first calculating 
                                        % derivatives of the mapping from
    r = sqrt(d2(:,i));                  % fi and th to x,y,z.
    du_dr = dfi_dr(r);
        
    du_dx = du_dr.*(xyz(i,1)-xyz(:,1))./r;
    du_dy = du_dr.*(xyz(i,2)-xyz(:,2))./r;
    du_dz = du_dr.*(xyz(i,3)-xyz(:,3))./r;
    du_dx(i) = 0; du_dy(i) = 0; du_dz(i) = 0;

    du_dfi = du_dx*dx_dfi + du_dy*dy_dfi + du_dz*dz_dfi;
    du_dth = du_dx*dx_dth + du_dy*dy_dth + du_dz*dz_dth;
        
    L_u = -(cos(al)-tan(xyz(i,5))*xyz(i,7)*sin(al))*du_dfi + ...
            sin(al)*xyz(i,6)*du_dth;    % L-operator at node i evaluated 
                                        % for the different RBFs
    D(i,:) = L_u'*la;                   % Row i of the DM computed
end

fig1 = figure(1)
plot(eig(D),'k.'); axis square          % Plot the eigenvalues of the DM
title('Eigenvalues of Differentiation Matrix: Global RBF')
% D: Differentiation Matrix needed.
D;


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
fig2 = figure(2)
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

% % % % -------------------------------------------------------------------------
% % % % Plots of Time-Revolution: JSurf (Matlab Central FileID: #8763)
% % % 
% % % fig5 = figure(5)
% % % [X,Y] = meshgrid(ptsX,ptsY); 
% % % % Compute Z
% % % temp = (X.^2 + Y.^2)<=1;
% % % Z = sqrt(ones(size(X)) - temp.*(X.^2+Y.^2));   % How to determine SIGN???
% % % cosTheta = sqrt(ones(size(Z)) - Z.^2);
% % % cosPhi = X./cosTheta;
% % % r = acos(cosTheta.*cosPhi);
% % % ratio = r/R;
% % % indx = (ratio<1);
% % % Hval  = h0/2*(1+cos(pi.*ratio)).*indx;        
% % % 
% % % 
% % % jsurf(X,Y,Z,Hval); 
% % % axis vis3d; 
% % % axis off; 
% % % set(gca,'DataAspectRatio',[1 1 1]);
