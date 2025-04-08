% This script vaildates the methods to find the stationary solution of Markov Chain

close all
clear all
clc

%.. PATH
mfilepath = pwd;
idcs = strfind(mfilepath,'\');
libdir = mfilepath(1:idcs(end));
addpath([libdir, 'CommonSource'])

state = 0; % 0: hot, 1: cold
iter_max = 1000; % # of iteration
Xset = zeros(iter_max,1); % Var. for saving

% Run Simulation
rnd_set = rand(iter_max,1); % Better not to call r.n.g every iteration if possible
for ii = 1:iter_max
    rnd = rnd_set(ii);
    if state == 0 % State = Hot
        if rnd < 0.2 % Hot -> Hot
            state = 0;
        else % Hot -> Cold
            state = 1;
        end
    elseif state == 1 % State = Cold
        if rnd < 0.6 % Cold -> Hot
            state = 0;
        else % Cold -> Cold
            state = 1;
        end
    end
    Xset(ii) = state;
end

% Count the number
cnt_cold = sum(Xset);
cnt_hot = iter_max - cnt_cold;
pi_sim = [cnt_hot; cnt_cold];
pi_sim = pi_sim/sum(pi_sim);

% Long term behavior: eigenvalue
P = [0.2, 0.8; 0.6, 0.4]'; % Transpose since pi is defined to be column vector
[V,D] = eig(P); % Find eigenvector having eigenvalue == 1
pi_anl1 = V(:,2)/sum(V(:,2));

% Long term behavior: state space reduction
P100 = P^100;
pi_anl2 = limitdist(P); % Row vector form is used

[pi_sim, pi_anl1, pi_anl2]

% 
time_step = 500;
pi0 = [1; 0];
pi_set = zeros(2,time_step);
pi_set(:,1) = pi0;
for i = 1:time_step-1
    pi_set(:,i+1) = P*pi_set(:,i);
end
figure(1); hold on; grid on
plot(1:time_step, pi_set(1,:),...
     1:time_step, pi_set(2,:))
xlabel('Time Step')
ylabel('State Prob')
legend('State 1', 'State 2')

%% Test Large Matrix
n = 500;
P = rand(n,n);
for i = 1:n
    P(:,i) = P(:,i)/sum(P(:,i));
end

iter_max = 10;
tic
for i = 1:iter_max
    x3 = limitdist(P,3); % State-space reduction method
end
T3 = toc

tic
for i = 1:iter_max
    x2 = limitdist(P,2); % Built-in linear solver 
end
T2 = toc

tic
for i = 1:iter_max
    x1 = limitdist(P,1); % Power Iteration
end
T1 = toc

norm(x2-x1)
norm(x2-x3)