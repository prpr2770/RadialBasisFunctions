% APPM 7440: HW#3
% Question 1: Double Precision Arithmetic

function feval = evalRBFinterpolation(Z_data, Z_grid,lambda, epsilon, RBFtype)
%{
Z_data      : Data point vectors
Z_grid      : Grid point vectors
lambda      : RBF Interpolation Coefficients
epsilon     : Shape Parameter
RBFtype     : Type of RBF used
%}

% MQ-RBF/ GA-RBF methodology for function interpolation
% and computing error.
% --------------------------
% Compute as a function of shape-parameter
% --------------------------
% meshgrid and ndgrid to compute the A matrix
% --------------------------

[r,c] = size(Z_data);
if c ==1 % 1-Dimensional Data
    X_data = Z_data(:,1)';
    X_grid = Z_grid(:,1)';
    
    G = ndgrid(X_grid, X_data);
    GridX = G';
    DataX = ndgrid(X_data, X_grid);
    
    dX = DataX - GridX;
    
    R2 = dX.^2;
    
elseif c == 2 % 2-dimensional dataset
    X_data = Z_data(:,1)'
    Y_data = Z_data(:,2)'

        
    X_grid = Z_grid(:,1)'
    Y_grid = Z_grid(:,2)'

   
    % ==================================
    
    G = ndgrid(X_grid, X_data);
    GridX = G'
    DataX = ndgrid(X_data, X_grid)
    
   
    G = ndgrid(Y_grid, Y_data);
    GridY = G'
    DataY = ndgrid(Y_data, Y_grid)
    
    % ============================
    % Compute RBF
    
    dX = DataX - GridX;
    dY = DataY - GridY;
    R2 = dX.^2 + dY.^2;
    
else
    fprintf('>2-Dimensions NOT supported!')
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
        A = 1./(ONE + epsilon^2*R2).^(1/2);
        
    case 'IQ' % RBF: Inverse Quadratic
        fprintf('Inverse Quadratic \n')
        A = 1./(ONE + epsilon^2*R2);
        
    otherwise % RBF: Multi-Quadrics
        fprintf('Multi-Quadrics \n')
        A = (ONE + epsilon^2*R2).^(1/2);
end

feval = A*lambda

return