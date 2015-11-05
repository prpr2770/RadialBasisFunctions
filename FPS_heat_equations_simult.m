function FPS_heat_equations_simult
% 1D-Heat Equation: Fourier Pseudospectral
% PDE u_t = u_xx
% IC u(x,0) = 0.05/(1.05-cos(pi*x))
% for time t=0 to t=0.2
clear all
close all
clc

t0 = 0;tf = 0.2;
n = 32;                                 % set space resolution

% At every point in space, we run an ODE in time. 
x = linspace(-1,1,n+1); x(end) = [];    %layout datapoints; skip last.
u = (0.05./(1.05-cos(pi*x)) )';        % give Init. cond.

u0 = (u + i*u)';        % Give I.C of two simult eqs as complex no.s

v = (-pi^2*[0:n/2, n/2-1:-1:1].^2)';    % vector to multiply in fourier space
                                        % to get second derivatives.

% use ode45 to advance in time
[t,u] = ode45(@dudt, t0:0.01:tf, u0, [], v, x);    

u(:,n+1) = u(:,1);                      % Fill in right-edge in graphics.

fig1 = figure(1);
% Visualize PDE in Real Domain
% fprintf('size u: %f %f \n', size(u))
% fprintf('size t: %f  %f \n', size(t))
% fprintf('size [x,1]: %f  %f \n', size([x,1]))
mesh([x,1],t,real(u)); colormap([0,0,0])      


fig2 = figure(2);
% Visualize PDE in Imag Domain
% fprintf('size u: %f %f \n', size(u))
% fprintf('size t: %f  %f \n', size(t))
% fprintf('size [x,1]: %f  %f \n', size([x,1]))
mesh([x,1],t,imag(u)); colormap([0,0,0])      



function uxx = dudt(t,u,v, x)
%{
Calculate uxx by Fourier Pseudospectral Method
The fourier-coefficients helps carry the spatial grid inter-dependency
between the ODE.
%}

uf = fft(u);
uf = uf.*v;
uxx = ifft(uf);


