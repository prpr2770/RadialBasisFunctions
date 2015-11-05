
function LPhi = Lfi(type,ep,r)

switch type
    case 'GA'
    % Gaussian RBF
    LPhi =  ( 4*ep^4*r.^2.*exp(-(ep*r).^2));      
    
    case 'MQ'
    % MQ
    LPhi =  -(ep*r).^2.*(1 + (ep*r).^2).^(-3/2);  
        
    case 'IQ'
    % IQ
    LPhi = (8*ep^4*r.^2.*(1 + (ep*r).^2).^(-3) ); 

end

end
