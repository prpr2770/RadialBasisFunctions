function Phi = fi(type,ep,r)

switch type
    case 'GA'
    % Gaussian RBF
    Phi = (exp(-(ep*r).^2));
    
    case 'MQ'
    % MQ
    Phi = (1 + (ep*r).^2).^(1/2);    
        
    case 'IQ'
    % IQ
    Phi = (1 + (ep*r).^2).^(-1);    

end

end


    
    
    
       
    
    
       