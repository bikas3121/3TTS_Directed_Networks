function [init] =  InitialConditions(nCk)
% The function InitialConditions generates initial conditions for the whole
% network. 
% The function normrnd has the following syntax:
% normrnd(a, b, [m,n])
% It generates the normally distributed random numbers with mean 'a' and
% the standard deivation 'b' around the mean with [m,n] matrix size.

%% Input
% nCk  = [3 4 5 6];

%% Initial Condition Generation
n = sum(nCk);
init = [];
mean = 50;
% a = [-50 50 -100 100]; % mean of the initial conditions for each cluster
b = 30; % standard deviation
% means = zeros(length(nCk),1);
for i = 1:length(nCk)
%     a = (-1)^i*mean*(1/);
    a = (-1)^(i)*mean; 
    mn = a*(1/2*i);
    rng ('default');
    init_i = normrnd(mn,b, [2*nCk(i),1]);
%     init_i = normrnd(a(i)*2,b, [2*nCk(i),1]);
    init = [init; init_i];
end
end