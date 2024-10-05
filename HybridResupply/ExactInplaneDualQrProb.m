% Def
% Compute Exact In-Plane State Distrubution using Different Method
%
% Input

%
% Output

%
% Seungyeop, Han 

function [x, x_aft, x_bf, x_dmd] = ExactInplaneQrProb(type, n, Pd, Pav, Q, R)
    %.. State
    x = Q+R - (0:1:n-1);
    
    switch type
        case 1 % 100% Parking Available
            [x_aft, x_bf, x_dmd] = PerfectAvailable(n, Pd, Q, R);
            
        case 2 % Re-Order with Pav Probability
            [x_aft, x_bf, x_dmd] = BinaryAvailable(n, Pd, Pav, Q, R);
            
        case 3 % Re-Order with Pav(0), Pav(1), ..., Pav(i), ... Probability
            [x_aft, x_bf, x_dmd] = GeneralAvailable(n, Pd, Pav, Q, R);
    end

end


function [x_aft, x_bf, x_dmd] = PerfectAvailable(n, Pd, Q, R)
    
    %.. Check State Dimension
    if mod(n,Q) ~= 0
        error('Dimension of n is the multiple of Q');
    end
    
    %..
    m = n/Q;
    
    %..
    U = zeros(Q, n);
    for i = 1:Q
        U(i,i:end) = Pd(1:n+1-i);
    end
    
    %.. Distribution After Resupply
    x_aft = ones(1, Q)/Q;
    
    % Distribution Before Resupply
    x_bf = x_aft*U;
    x_aft = [x_aft, zeros(1,n-Q)];
    
    % Demand Distribution
    X_bf = reshape(x_bf,Q,m);
    x_dmd = sum(X_bf);
end

function [x_aft, x_bf, x_dmd] = BinaryAvailable(n, Pd, Pav, Q, R)

    %.. Check State Dimension
    if mod(n,Q) ~= 0
        error('Dimension of n is the multiple of Q');
    end
    
    %
    m = n/Q;
    k1 = Pav;
    k2 = 1-Pav;
    
    % U
    U = zeros(n, n);
    for i = 1:n
        U(i,i:end) = Pd(1:n+1-i);
    end
    
    % U0
    U0 = U(1:Q,1:Q);
    
    % B1, B2, B3...
    B = zeros(Q,Q,m-1);
    for i = 1:m-1
        idx = Q*i + 1;
        B(:,:,i) = U(1:Q,idx:idx+Q-1);
    end
    
    % Distribution After Resupply
    X_aft = zeros(m,Q);
    y = ones(1,Q)/Q;
    X_aft(1,:) = (k1*y)/(eye(Q)-k2*U0); % k1*y*inv(eye(Q)-k2*U0)
    for i = 2:m
        tmp = zeros(1,Q);
        for j = 1:i-1
            tmp = tmp + X_aft(j,:)*B(:,:,i-j);
        end

        X_aft(i,:) = (k2*tmp)/(eye(Q)-k2*U0);
    end
    x_aft = X_aft';
    x_aft = (x_aft(:))';

    % Distribution Before Resupply
    x_bf = x_aft*U;
    
    % Demand Distribution
    X_bf = reshape(x_bf,Q,m);
    x_dmd = sum(X_bf);
end

function [x_aft, x_bf, x_dmd] = GeneralAvailable(n, Pd, Pavd, Q, R)

    m = n/Q;
    
    % P_plane-
    Ptrs = zeros(n, n);
    for i = 1:n
        Ptrs(i,i:end) = Pd(1:n+1-i);
    end
    
    % P_plane+
    Pre = zeros(n, n);
    for i = 1:m
        for j = 1:i
            x_idx = (Q*(i-1) + 1):(Q*i);
            y_idx = (Q*(j-1) + 1):(Q*j);
            
            %.. Prob
            if j == 1
                p_ij = Pavd(i);
            else
                p_ij = 1;
                for k = 1:j-1
                    p_ij = p_ij*(1-Pavd(i-k+1));
                end
                p_ij = p_ij*Pavd(i-j+1);
            end
            
            Pre(x_idx,y_idx) = p_ij*eye(Q);
        end
    end
    
    %.. P Plane and stationary dist
    Pp = Pre*Ptrs;
    for i = 1:n
        Pp(i,:) = Pp(i,:)/sum(Pp(i,:));
    end
    x_bf = limitdist(Pp);
    x_aft =  x_bf*Pre;
    

    % Demand Distribution
    X_bf = reshape(x_bf,Q,m);
    x_dmd = sum(X_bf);
end


function p=limitdist(P)
%Obtain the stationary probability distribution
%vector p of an irreducible, recurrent Markov
%chain by state reduction. P is the transition
%probabilities matrix of a discrete-time Markov
%chain or the generator matrix Q.
% https://www.math.wustl.edu/~feres/Math450Lect04.pdf

[ns ms]=size(P);
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
p = p';
end
