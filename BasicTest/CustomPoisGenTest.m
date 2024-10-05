% This script tests the custom-made poisson random number generator
% To validate the custom-fucntion, you will need the toolbox for poissrnd function

clear all
close all
clc

%.. PATH
mfilepath = pwd;
idcs = strfind(mfilepath,'\');
libdir = mfilepath(1:idcs(end));
addpath([libdir, 'CommonSource'])

%.. Test Param
lamda = 2; % Dist. Parameter
n = 1000; % Number of function call

tic
for i = 1:n
    Zmy = CustomPoisRnd(lamda,10000,1);
end
Tmy = toc;

tic
for i = 1:n
    Zmat = poissrnd(lamda,10000,1);
end
Tmat = toc;

disp([Tmy, Tmat]/n)

edges = 0:1:11;
h1 = histcounts(Zmy,edges);
h2 = histcounts(Zmat,edges);

figure(1)
bar(edges(1:end-1),[h1; h2])