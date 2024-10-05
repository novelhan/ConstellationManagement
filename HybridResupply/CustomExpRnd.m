function [Y] = CustomExpRnd(mu, m, n)
    % Input Handling
    if nargin == 1
        m = 1;
        n = 1;
    elseif nargin == 2
        n = 1;
    end
    
    % Mean = 1/lambda
    p = 1/mu;
    
    % Generate uniform distribution
    Z = rand(m,n);
    
    % Apply Inverse Transform
    Y = -1/p*log(1-Z);
end