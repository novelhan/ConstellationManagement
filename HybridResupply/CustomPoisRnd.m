function [Y] = CustomPoisRnd(p, m, n)
    % Input Handling
    if nargin == 1
        m = 1;
        n = 1;
    elseif nargin == 2
        n = 1;
    end
    
    % Compute Poisson CDF (consider upto mean + 10-sigma)
    Zset = 0:1:ceil(11*p);
    PDF = p.^Zset*exp(-p)./factorial(Zset); 
    CDF = cumsum(PDF);
    CDF(end) = 1;
    
    % Generate uniform distribution
    Z = rand(m,n);
    Y = zeros(m,n);
    for i = 1:m
        for j = 1:n
            tmp = CDF - Z(i,j);
            idx = find(tmp>=0,1);
            Y(i,j) = Zset(idx);
        end
    end
end