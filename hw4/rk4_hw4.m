function [t,u] = rk4_hw4(fun_dudt_hw4, domain, initCond, spatialGridLayout )
'''Script implementing rk4 algo to solve ODE over t'''

% time-series over which to compute ODE
t = domain; 
dt = t(2)- t(1);        % Time-Grid Spacing

% space-vector at all spatial points
u0 = initCond;                  % columnVector obtained at DifferentSpatial points
u = zeros(length(u0),length(t));
u(:,1) = u0;

for i = 2:numel(t)
    T = t(i-1);
    U = u(:,i-1);
    k1 = dt*fun_dudt_hw4(T,U, spatialGridLayout);
    k2 = dt*fun_dudt_hw4(T+dt/2, U + k1/2, spatialGridLayout);
    k3 = dt*fun_dudt_hw4(T+dt/2, U + k2/2, spatialGridLayout);
    k4 = dt*fun_dudt_hw4(T+dt, U + k3, spatialGridLayout);
    % update value of function at next-time step
    u(:,i) = U + 1/6*(k1+2*k2+2*k3+k4);
end
