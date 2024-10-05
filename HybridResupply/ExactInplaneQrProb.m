% Def
% Compute Exact In-Plane State Distrubution using Different Method
%
% Input

%
% Output

%
% Seungyeop, Han 

function [pi_p, Pav, Pdav] = ExactInplaneQrProb(f, Q, R, mu, m)
    
    %.. State
    xmax = Q + R;
    x = 0:1:xmax;
    nx = length(x);
    
    %.. Failure Transition Matrix
    Pf = zeros(nx);
    for i = 1:nx
        Pf(i:end,i) = poisspdf(1:nx+1-i, nx*f)';        % Constant Failure Rate
%         Pf(i:end,i) = poisspdf(1:nx+1-i, (nx+1-i)*f)';  % Stock Level Dependant Failure Rate
        Pf(end,i) = 1 - sum(Pf(1:end-1, i));
    end
    U0 = Pf(1:Q, 1:Q);
    U1 = Pf(Q+1:end, Q+1:end);
    B = Pf(Q+1:end, 1:Q);
    
    %.. Resupply Transition Matrix
    Pr = [ [eye(Q); zeros(R+1, Q)], [eye(R+1); zeros(Q,R+1)] ];
    
    
    if R+1 == Q
        L = eye(R+1);
    else
        L = [eye(R+1), zeros(R+1, Q-R)];
    end
    
    %.. Stationary Matrox (May optimize the code)
    C0 = (eye(Q) - U0)\B;
    C1 = (1 - exp(-mu))*U1^m*((eye(R+1) - exp(-mu)*U1)\L);
    Plim = C0*C1;
    
    C2 = eye(R+1);
    for i = 1:m-1
        C2 = C2 + U1^i;
    end
    C2 = C2 + U1^m/(eye(R+1) - exp(-mu)*U1);
    
    %.. Stationary Distribution
    pi_s = limitdist(Plim);
    pi_r = pi_s*C0;
    pi_0 = pi_s/(eye(Q) - U0);
    pi_0 = pi_0/sum(pi_0);
    pi_1 = pi_r*C2;
    pi_1 = pi_1/sum(pi_1);
    
    %.. Full Stationary Distribution
    avgPd = x*Pf';
    avgPir = (Q:-1:1)*pi_s';
    N0 = 0.5 + avgPir/avgPd;
    N1 = 1 + m + 1/(exp(mu)-1);
    b0 = N0/(N0+N1);
    b1 = N1/(N0+N1);
    pi_p = [b0*pi_0, b1*pi_1];
    
    %.. Parking Availability
    Pav = 0;
    for i = 2:nx-1
        Pav = Pav + Pf(i)*sum(pi_p(1:Q+R+2-i));
    end
    Pav = Pav/(1-Pf(1));
    
    Pdav = zeros(1,nx);
    for i = 1:nx
        Pdav(i) = sum(pi_p(1:Q+R+2-i));
    end
    
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