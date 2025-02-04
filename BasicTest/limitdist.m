function p = limitdist(P,method)

%.. Select computation method
if nargin == 1
    method = 1;
end
    
if method == 2
    % State Space Reduction Method
    % Obtain the stationary probability distribution vector p of an irreducible, recurrent Markov
    % chain by state reduction. P is the transition probabilities matrix of a discrete-time Markov
    % chain or the generator matrix Q.
    % https://www.math.wustl.edu/~feres/Math450Lect04.pdf
    
    [ns, ~]=size(P);
    n=ns;
    p=zeros(n);
    while n>1
        n1=n-1;
        s=sum(P(n,1:n1));
        P(1:n1,n)=P(1:n1,n)/s;
        n2=n1;
        while n2>0
            P(1:n1,n2)=P(1:n1,n2)+P(1:n1,n)*P(n,n2);
            n2=n2-1;
        end
        n=n-1;
    end
    %backtracking
    p(1)=1;
    j=2;
    while j<=ns
        j1=j-1;
        p(j)=sum(p(1:j1).*(P(1:j1,j))');
        j=j+1;
    end
    p=p/(sum(p));
    
else
    % Simple solution approach for linear equation
    % Rmk: Based on the several test, this approach has better computational speed w/ same accuracy
    [n, ~] = size(P);
    A = [eye(n)-P; ones(1,n)];
    b = [zeros(n,1); 1];
    p = A\b;
end

end