% Assignment 5: 
% Prasanth Prahladan(100817764)
% Submitted: 4 Nov, 2014
% ==============================

%P.1. 

A = [2.7 -15.1 -1.2; 0 0 1; 1 0 0]
x0 = [0.2; -0.1; 1]


T1 = expm(A)
T2 = expm(2*A)
T3 = expm(3*A)


x1 = T1*x0
x2 = T2*x0
x3 = T3*x0

eig(A)

%========================================================
%P.2
%program to print phase trajectory of given 2-var nonlinear system
%function phaseTrajectory(X0,Y0)
diffeq =@(t,X)[X(2); (1-X(1).*X(1)).*X(2) - X(1)];
x1 = linspace(-3,3,100);
x2 = linspace(-3,3,100);
[x,y] = meshgrid(x1,x2);
u = zeros(size(x));
v = zeros(size(y));
t=0;
for i = 1:numel(x)
   
    Xprime = diffeq(t,[x(i);y(i)]);
    u(i) = Xprime(1);
    v(i) = Xprime(2);
end

quiver(x,y,u,v,'r'); figure(gcf)


%print trajectories for inital points at (x1,x2) at 
% random points starting in (x,y) \in [-1,1]x[-1,1]
figure(1)
x1_orig = 2*rand(1,2)+ [-1 -1]
x2_orig = 2*rand(1,2)+ [-1 -1]
hold on
for x10 = x1_orig
    hold on
    for x20 = x2_orig
        [ts,xs] = ode45(diffeq,[0,100],[x20;0]);
        plot (xs(:,1),xs(:,2),'g')
        plot (xs(1,1),xs(1,2),'ro')
        plot (xs(end,1),xs(end,2),'ks')
        [ts,ys] = ode45(diffeq,[0,100],[0;x20]);
        plot (ys(:,1),ys(:,2),'b')
        plot (ys(1,1),ys(1,2),'ro')
        plot (ys(end,1),ys(end,2),'ks')
    end
    hold off
end

% In the plot, the points with CIRCLE mark starting points. 
% The trajectories terminate at the SQUARE mark. 
% We observe that the Equilibrium point (0,0) is UNSTABLE, as the
% trajectories move away from it onto the Limit Cycle.

%========================================================
% Q.3.A)
A=[1 -1;1 -2];
B=[1;0.2];
K = place(A,B,[-1,-2])
%========================================================
% Q.3.B)
% Simulink model.