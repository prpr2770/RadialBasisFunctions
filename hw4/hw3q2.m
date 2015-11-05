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

%%
% Choose the data-set you want from ABOVE
% (x2,y2): Haltonset ; (x1,y1): RandomNumbers

% #COLUMNS = # Dimensions of the Spatial Vectors
% #ROWS = # Data Points

x0 = sym(x2);
y0 = sym(y2);

% Divide Data Points into GRID(TRAIN) and TEST LOCATIONS
% Grid-Training Data
x_trn = x0(1:N/2)';
y_trn = y0(1:N/2)';
Z_trn = [x_trn; y_trn]';
Z_trn = sym(Z_trn);

% Testing Data
x_tst = x0(N/2+1:end)';
y_tst = y0(N/2+1:end)';
Z_tst = [x_tst; y_tst]';
Z_tst = sym(Z_tst);

%%  MAIN SECTION: Iterating over SHAPE PARAMTER

% Describe Test Function
testFunc = @(Z)( sym(59) ./ sym(67 + ( sym(Z(:,1)) + 1/sym(7) ).^2 + ( sym(Z(:,2)) - 1/sym(11) ).^2 ) );

% Evaluate Test Function : Train Data
fvals_trn = testFunc(Z_trn);
fvals_trn = sym(fvals_trn);
% fprintf('\nSize fvals_trn \t'); size(fvals_trn)

% Evaluate Test Function : Testing Data
fvals_tst = testFunc(Z_tst);
fvals_tst = sym(fvals_tst);
% fprintf('\nSize fvals_tst \t'); size(fvals_tst)

% ===============================================
% Evlauation for Training Data

X_trn = sym(Z_trn(:,1));
Y_trn = sym(Z_trn(:,2));

% grid ------------------------
[gx, gy] = meshgrid(X_trn,Y_trn);
GY = sym(gy');
GX = sym(gx);

% Compute RBF
dX = sym(GX' - GX);
dY = sym(GY' - GY);
R2Trn = sym(dX.^2 + dY.^2);
rTrn = sym(sqrt(R2Trn));
% fprintf('\nSize rTrn \t');size(rTrn)

% ========================================
% Computations for Testing Data
X_tst = sym(Z_tst(:,1));
Y_tst = sym(Z_tst(:,2));

[nodeX,gridX] = ndgrid(X_tst,X_trn);
[nodeY,gridY] = ndgrid(Y_tst,Y_trn);

R2tst = (sym(nodeX)-sym(gridX)).^2 + (sym(nodeY) - sym(gridY)).^2;
rTst = sqrt(sym(R2tst));
% fprintf('\nSize rTst \t');size(rTst)
% ===============================================
% ------------------------------------------------------------
% Epsilon = ones(size(Epsilon)); % Test Case for Accuracy
% When eps = 1, the Reciprocal Condition-Number should be closer to 1.
% ------------------------------------------------------------
Epsilon = 0.001:0.001:1;

error = zeros(size(Epsilon));
Condition = zeros(size(Epsilon));
Lambda = zeros(size(Epsilon));
Fvals = zeros(size(Epsilon));
% ===============================================
% Description of RBF Function
fi = @(r,ep)(exp(-(ep*r).^2)); % GA-RBF

ACC = 2^10; % VPA Precision Digits
profile ON
% ----------------------------------------------------------------
for count = 1:length(Epsilon)
% count = 1;
    
    epsilon = Epsilon(count);
    fprintf('epsilon \t %f \n',epsilon)
    epsilon = sym(epsilon);
    
    % Compute Lambda
    % --------------------------
    A_trn = fi(rTrn,epsilon);
    
    nmA_trn = vpa(A_trn,ACC);
%     fprintf('\nSize A_trn \t');size(A_trn)
    
%     invA_trn = vpa(inv(A_trn),ACC);
%     fprintf('\nSize invA_trn \t');size(invA_trn)    
    
%     fprintf('\nSize fvals_trn \t');size(fvals_trn)
    
    lambda = nmA_trn\fvals_trn;
    lambda = sym(lambda);
%     fprintf('\n Size lambda \t');size(lambda)
    
    % Evaluate the function at Test Points
    % --------------------------
    A_tst = fi(rTst,epsilon);
%     A_tst = sym(A_tst);
%     fprintf('\nSize A_tst \t');size(A_tst)
    
    feval_tst = A_tst*lambda;
    feval_tst = sym(feval_tst);
%     fprintf('\nSize feval_tst \t');size(feval_tst)
    
    % Compute the parameters to plot.
    % --------------------------
    Condition(count) = cond(nmA_trn);
    err = vpa(feval_tst - fvals_tst, ACC);
    error(count) = max(abs(err));
    
    % unnecessary params
    Lambda(count) = max(abs(vpa(lambda,ACC)));
    Fvals(count) = max(abs(vpa(fvals_trn,ACC)));
    
end

profile VIEWER

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