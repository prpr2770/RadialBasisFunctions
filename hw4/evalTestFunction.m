% APPM 7440: HW#3
% Question 1: Double Precision Arithmetic

function fvals = evalTestFunction(Z)
x = Z(:,1); % Row vector : X-dim
y = Z(:,2); %
fvals = 59 ./ (67 + (x + 1/7).^2 + (y - 1/11).^2);

