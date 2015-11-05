% APPM 7440: HW#3
% Question 1: Double Precision Arithmetic
% RBF - DIRECT 
% ========================================================================
close all
clear all
clc

N = 42;         % Total sample points needed : Always Ensure EVEN!!

%%
% GENERATE POINTS WITHIN UNIT CIRCLE - can be done as a function.

% plot unit circle
theta = 0:0.01:2*pi;
circx = cos(theta);
circy = sin(theta);
fig1 = figure(1)
plot(circx,circy,'r')
hold on
%%
% Generate Uniformly Random distributed points in Unit Disc

rndTheta = 2*pi*rand(N,1);
rndRad1 = sqrt(rand(N,1));
x1 = rndRad1.*cos(rndTheta);
y1 = rndRad1.*sin(rndTheta);
scatter(x1,y1)
title('Uniformly Distributed Points in Unit Disc')
hold off

% xvecs = [x1' y1']; % Nx2 matrix.

%%
% Generate Random points using HaltonSet

rng default
p = haltonset(2,'Skip',1e3,'Leap', 1e2);
p = scramble(p,'RR2');      % Reverse radix scrambling
X0 = net(p,N);
X0 = -1*ones(size(X0)) + 2*X0;
R = sqrt(X0(:,1).^2 + X0(:,2).^2);
inds = find((R<=1));
x2 = X0(inds,1);
y2 = X0(inds,2);

% -----------------------------------
fig2 = figure(2)
% subplot(2,1,1)
scatter(x2,y2,'b')
title('Haltonset Points in Circle')
hold on
plot(circx,circy,'r')
hold off
% subplot(2,1,2)
% Rfin = R.*(R<=1);
% scatter(Rfin, Rfin)
% -----------------------------------


% How to ensure Even Number of Data Points
% --- DO THIS LATER --



%%
% Choose the data-set you want from ABOVE
% (x2,y2): Haltonset ; (x1,y1): RandomNumbers

% #COLUMNS = # Dimensions of the Spatial Vectors
% #ROWS = # Data Points

x0 = x2;
y0 = y2;

% x = floor(10*rand(1,5));
% y = floor(10*rand(1,5));

% Divide Data Points into GRID(TRAIN) and TEST LOCATIONS
x_trn = x0(1:N/2)';
y_trn = y0(1:N/2)';
Z_trn = [x_trn; y_trn]';   
size(Z_trn)

x_tst = x0(N/2+1:end)';
y_tst = y0(N/2+1:end)';
Z_tst = [x_tst; y_tst]';


%%
% Checking for errors when the Test Data is same as Training Data
% Z_tst = Z_trn(2:4,:);
 Z_tst = Z_tst;
%%  MAIN SECTION: Iterating over SHAPE PARAMTER

% Evaluate Test Function : Train Data
fvals_trn = evalTestFunction(Z_trn);
size(fvals_trn)

% Evaluate Test Function : Testing Data
fvals_tst = evalTestFunction(Z_tst);
size(fvals_tst)

Epsilon = 0.01:0.001:2;
% ------------------------------------------------------------
% Epsilon = ones(size(Epsilon)); % Test Case for Accuracy
% When eps = 1, the Reciprocal Condition-Number should be closer to 1.
% ------------------------------------------------------------
error = zeros(size(Epsilon));
Condition = zeros(size(Epsilon));
Lambda = zeros(size(Epsilon));
Fvals = zeros(size(Epsilon));

% ----------------------------------------------------------------
for count = 1:length(Epsilon)
% Compute Lambda
% --------------------------
% fprintf('Z_trn Dimensions \n')
% size(Z_trn)
epsilon = Epsilon(count);
A_trn = getRBFmatrix(Z_trn,epsilon,'GA');
% fprintf('Size A_trn \n')
% size(A_trn)
lambda = A_trn\fvals_trn;

% Evaluate the function at 
feval_tst = evalRBFinterpolation(Z_tst,Z_trn,lambda, epsilon,'GA');

fprintf('Size feval_tst \t')
size(feval_tst)
fprintf('Size fvals_tst \t')
size(fvals_tst)

% Compute the parameters to plot.
Condition(count) = cond(A_trn);
error(count) = max(abs(feval_tst - fvals_tst));

% unnecessary params
Lambda(count) = max(abs(lambda));
Fvals(count) = max(abs(fvals_trn));

% break         % Used for Testing!

end

count
fig4 = figure(4)
subplot(2,1,1)
% plot(Epsilon, Lambda)
plot(log10(Lambda))        % TEST
xlabel('Epsilon')
ylabel('Log Max Lambda')
subplot(2,1,2)
% plot(Epsilon, Fvals)
plot(Fvals)         % TEST
xlabel('Epsilon')
ylabel('Fvals')

fig5 = figure(5)
subplot(2,1,1)
plot(Epsilon, log10(error))
% plot(error)         % TEST
title('Error Accuracy')
xlabel('Epsilon')
ylabel('Log Error')
subplot(2,1,2)
plot(Epsilon, log10(Condition))
% plot(Condition)     % TEST
title('Condiiton Number of A-Matrix')
xlabel('Epsilon')
ylabel('Log cond')