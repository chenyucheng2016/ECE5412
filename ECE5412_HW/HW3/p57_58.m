%p57,p58
clc,clear,close all
%generate observations Y
a = 5;
b = 0.5;
omega_v = .5;
omega_n = 1;
x0 = normrnd(0,1,1);
T = 40;
Y = zeros(1,T);
x_cur = x0;
for i = 1:T
    x_cur = x_cur + a - b*exp(x_cur) + normrnd(0,omega_v); 
    y = x_cur + normrnd(0,omega_n);
    Y(i) = y;
end

%estimate the parameters a
a_prior = normrnd(2,2);
[X,x0] = sampleX(a_prior,b,T,omega_v);
N = 1000; %gibbs sampling iterations
%Gibbs sampling loop
As = zeros(1,N);
Bs = zeros(1,N);
for i = 1:N
    paramSamples = MCMC_params(X,Y,omega_v,omega_n,x0);
    param = paramSamples(:,unidrnd(size(paramSamples,2)))
    [X,x0] = SMCTraj(param(1),b,omega_v,omega_n,Y,T)
    As(i) = param(1);
    Bs(i) = param(2);
end
hist(As,200);
xticks(-15:1:10)
function [X,x0] = SMCTraj(a,b,omega_v,omega_n,Y,T)
N = 8000;
x0 = normrnd(0,1,[N,1]);
sampleTrajs = zeros(N,T+2);
sampleTrajs(:,2) = x0;
for j = 1:N
    x_cur = x0(j);
    for i = 1:T
        x_cur = x_cur + a - b*exp(x_cur) + normrnd(0,omega_v);
        sampleTrajs(j,i+2) = x_cur;
        y = Y(i);
        sampleTrajs(j,1) = sampleTrajs(j,1) + log(max(normpdf(y,x_cur,omega_n),exp(-100)));
    end

end
    [~,sampleTrajIndex] = max(sampleTrajs(:,1));
    x0 = sampleTrajs(sampleTrajIndex,2);
    X = sampleTrajs(sampleTrajIndex,3:end);
end
function [X,x0] = sampleX(a,b,T,omega_v)
X = zeros(1,T);
x0 = normrnd(0,1,1);
x_cur = x0;
for i = 1:T
    x_cur = x_cur + a - b*exp(x_cur) + normrnd(0,omega_v); 
    X(i) = x_cur;
end
end

function params = MCMC_params(X,Y,omega_v,omega_n,x0)
ab0 = normrnd(3,2);
N = 8000;
probVar = 0.01;
%a = a0;b = b0;
a_prev = ab0(1);b_prev = 0.5;%b_prev = ab0(2);
params = zeros(2,N);
%params(:,1) = [a0;b0];
for i = 1:N
    %ab = mvnrnd([a_prev,b_prev],[probpVar,0;0,probpVar]);%proposal distribution
    a = normrnd(a_prev,probVar);
    %b = ab(2);
    %%%%
    b = 0.5;
    %%%
    logProb1 = evalTrajectoryLogProb(a_prev,b_prev,omega_v,omega_n,X,Y,x0);
    logProb2 = evalTrajectoryLogProb(a,b,omega_v,omega_n,X,Y,x0);
    logProb1 = logProb1 + normpdf(a,a_prev,probVar);
    logProb2 = logProb2 + normpdf(a_prev,a,probVar);
    if logProb2 > logProb1
        params(:,i) = [a,b];
        a_prev = a;
        b_prev = b;
    else
        u = rand();
        if u <= exp(logProb2-logProb1)
            params(:,i) = [a;b];
            a_prev = a;
            b_prev = b;
        else
            params(:,i) = [a_prev;b_prev];
        end
    end
end
end


function logProb = evalTrajectoryLogProb(a,b,omega_v,omega_n,X,Y,x0)
logProb = log(normpdf(x0,3,1));
x_prev = x0;
for i = 1:length(X)
    x = X(i);
    y = Y(i);
    logProb = logProb + log(max(normpdf(x,x_prev + a - b*exp(x_prev),omega_v),exp(-100))) + log(max(normpdf(y,x,omega_n),exp(-100)));
end
end