% APPM 7440: HW#3
% Question 1: Double Precision Arithmetic


function A = getRBFmatrix(Z,epsilon, RBFtype)
%{
Z           : Data point vectors
epsilon     : Shape Parameter
RBFtype     : Type of RBF used
%}

% 
% fprintf('Z - Grid Locations')
% Z

[rows,cols] = size(Z);
if cols ==1
    X = Z(:,1)'; % row vector
    [Ggrid, Gvar] = meshgrid(X,X);
    dX = Gvar - Ggrid;
    R2 = dX.^2;
    
elseif cols == 2 % 2-dimensional dataset
    % row-vectors
    X = Z(:,1)';
    Y = Z(:,2)';
    
    % grid ------------------------
    [gx, gy] = meshgrid(X,Y);
    GY = gy';
    GX = gx;
    
    % Compute RBF
    dX = GX' - GX;
    dY = GY' - GY;
    R2 = dX.^2 + dY.^2
    
else
    prinft('ERROR: Does not support greater than 2 Dimensions!')
end



% Create RBF Matrix
% ---------------------------
ONE = ones(size(R2));
switch RBFtype
    case 'GA' % RBF: Gaussian
        fprintf('Gaussian \n')
        A = exp(-(epsilon^2)*R2)
        
    case 'IMQ' % RBF: Inverse Multiquadric
        fprintf('Inverse Multiquadric \n')
        A = 1./(ONE + (epsilon^2)*R2).^(1/2);
        
    case 'IQ' % RBF: Inverse Quadratic
        fprintf('Inverse Quadratic \n')
        A = 1./(ONE + (epsilon^2)*R2);
        
    otherwise % RBF: Multi-Quadrics
        fprintf('Multi-Quadrics \n')
        A = (ONE + (epsilon^2)*R2).^(1/2);
        
end
% fprintf('A - matrix VALUE \n')
% A
%
% fprintf('inv(A) - matrix VALUE \n')
% inv(A)

return