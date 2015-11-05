function OneDim_HeatEquation
% 1D-Heat Equation: Fourier Pseudospectral
% RBF Course: HW#2
% Demo code of ode45/Fourier-PS solution of 
% PDE u_t = u_xx
% IC u(x,0) = 0.05/(1.05-cos(pi*x))
% for time t=0 to t=0.2
clear all
close all
clc


t0 = 0;
tf = 0.2;

n = 32;                                 % set space resolution

% At every point in space, we run an ODE in time. 
x = linspace(-1,1,n+1);     % layout datapoints
% x(end) = [];                % skip last
u0 = (0.05./(1.05-cos(pi*x)) )';        % give Init. cond.

% use rk4 to advance in time
[t,u] = rk4(@fun_dudt, t0:0.01:tf, u0, x);    

fig1 = figure(1)
mesh(x,t,u'); colormap([0,0,0])      %Plot result

function [t,u] = rk4(fun_dudt, domain, initCond, spatialGridLayout )
'''Script implementing rk4 algo to solve ODE over t'''

% time-series over which to compute ODE
t = domain; 
dt = t(2)- t(1);        % Time-Grid Spacing

% space-vector at all spatial points
u0 = initCond;                  % columnVector obtained at DifferentSpatial points
u = zeros(length(u0),length(t));
u(:,1) = u0;

for i = 2:numel(t)
    T = t(i-1)
    U = u(:,i-1);
    k1 = dt*fun_dudt(T,U, spatialGridLayout);
    k2 = dt*fun_dudt(T+dt/2, U + k1/2, spatialGridLayout);
    k3 = dt*fun_dudt(T+dt/2, U + k2/2, spatialGridLayout);
    k4 = dt*fun_dudt(T+dt, U + k3, spatialGridLayout)
    % update value of function at next-time step
    u(:,i) = U + 1/6*(k1+2*k2+2*k3+k4);
end


%{
function uxx = fun_dudt(t, U, X)
''' 1D-DiffEq ODE Function describing the dudt Velocity Vector '''
%{
U: Vector that contains values for U(x0:k:xF)
t:
X: Spatial Grid used
%}



X;                                      % Spatial-Grid 
dX = X(2) - X(1);                       % Spatial-Grid Spacing

Uxx = U;                                % Initialize Output    

% Determine the FD-Stencil
FD4_stencil = (1/dX^2)*[-1/12 4/3 -5/2 4/3 -1/12];
lenStencil = length(FD4_stencil);
halfStencil = floor((1/2)*length(FD4_stencil));

% % % % Extend Stencil for length of Data Vector
% % % ext_Stencil = zeros(1,length(U));       % Extended Stencil
% % % ext_Stencil(1:5) = FD4_stencil;
% % % ext_Stencil = circshift(ext_Stencil,[0 -halfStencil]);

shifted_U = circshift(U,[0 halfStencil]);

% Approximate the Second-Derivative using Spatial Grid Points
% The LIMITS of the LOOP is IMPORTANT!
for i = 1:numel(X)
    % how to handle the Boundary conditions????
%     fprintf('size : shifted_U \t %f %f \n',size(shifted_U));
if i< numel(X) - halfStencil -1 %0.5*length(FD4_stencil)
    u = shifted_U(i:i+4) ;
else
    
    diff = numel(X)-i
    u = zeros(lenStencil,1);
    length(shifted_U(i:end))
    length(U(1:lenStencil - diff -1))
    u(1:diff+1) = shifted_U(i:end)
    u(diff+2:end) = U(1:lenStencil - diff -1)
end
    fprintf('\n \n');
    fprintf('size : i \t%f\n',i);
    fprintf('size : u \t%f %f \n',size(u));
    fprintf('size : FD4_stencil \t%f %f \n',size(FD4_stencil));
    duxx = FD4_stencil*u;                      % Inner Product
    Uxx(i) = duxx;
    shifted_U = circshift(shifted_U, [0 -1]);
end

% Output final value
uxx = Uxx;
%}

function uxx = fun_dudt(t, U, X)
seed = 5;
rng(seed);
A = rand(length(U));
A = 0.5*(A + A');

uxx = A*U;

uxx = A*U;