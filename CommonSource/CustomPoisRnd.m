% Def
% Generate Poisson Random Number using 'rand' function
%
% Input
% p: Poisson Parameter
% m: Row size of the output
% n: Column size of the output

% Output
% Y: Generate Random Number

function [Y] = CustomPoisRnd(p, m, n)
    %.. To boost the computation speed, we precompute the factorial
    persistent FACTORIAL
    if isempty(FACTORIAL)
        FACTORIAL = [  1
                       1
                       2
                       6
                      24
                     120
                     720
                    5040
                   40320
                  362880
                 3628800
                39916800
               479001600
              6227020800
             87178291200];
    end

    % Input Handling
    if nargin == 1
        m = 1;
        n = 1;
    elseif nargin == 2
        n = 1;
    end
    
    % Compute Poisson CDF (consider upto mean + 10-sigma)
    Zset = 0:1:ceil(11*p);
%     PDF = p.^Zset*exp(-p)./factorial(Zset); 
    PDF = p.^Zset*exp(-p)./FACTORIAL(Zset+1); 
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