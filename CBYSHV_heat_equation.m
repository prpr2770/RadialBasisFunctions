function CBYSHV_heat_equation
% 1D-Heat Equation: Chebyshev-PS
% RBF Course: HW#2
% Demo code of ode45/Chebyshev-PS solution of 
% PDE u_t = u_xx
% IC u(x,0) = 0.05/(1.05-cos(pi*x))
% for time t=0 to t=0.2
clear all
close all
clc


t0 = 0;
tf = 0.2;

n = 32;                                 % set space resolution

% obtain 2nd-Derivative Matrix
g = [1 0 0; 1 0 0];
[x, D2t] = cheb2bc(n+1,g);
fprintf('Size x: %f %f\n', size(x))
CDM = D2t;
fprintf('Size CDM: %f %f\n', size(CDM))


% % Layout the Chebyshev Points
u0 = (0.05./(1.05-cos(pi*x)))';        % give Init. cond.
u0(end) = u0(1);
fprintf('Size u0: %f %f\n', size(u0))

fig1 = figure(1);
plot(x,u0)      % determine how function varies 

% use ode45 to advance in time
[t,u] = ode45(@dudt, t0:0.001:tf, u0, [], CDM, x);    


% Verify the following
% u(:,n+1) = u(:,1);                      % Fill in right-edge in graphics.

fig2 = figure(2);
fprintf('Figure BEGINS ================================= \n')
fprintf('size u: %f %f \n', size(u))
fprintf('size t: %f  %f \n', size(t))
fprintf('size x: %f  %f \n', size(x))
mesh(x',t,u); colormap([0,0,0])      %Plot result

function uxx = dudt(t,u,DM, x)
%{
Calculate uxx by Chebyshev Pseudospectral Method
%}
uxx = DM*u;


