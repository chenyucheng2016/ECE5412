clc,clear,close all
%question 2c
numSamples = 200000;
propVar = 5;
x = 1;
x_arr = zeros(1,numSamples);
for i = 1:numSamples
    x_cand = normrnd(x,propVar);
    x_prob = cos(x)^2*sin(2*x)^2*normpdf(x);
    x_cand_prob = cos(x_cand)^2*sin(2*x_cand)^2*normpdf(x_cand);
    if x_cand_prob > x_prob
        x = x_cand;
    else
        u = rand;
        if u < x_cand_prob/x_prob
            x = x_cand;
        end
    end
    x_arr(i) = x;
end

%state is z, 
%generate state and observation data
z = 1;
tspan = 500;
z_arr = zeros(1,tspan);
y_arr = zeros(1,tspan);
for i = 1:tspan
    z = z + normrnd(0,1);
    y = z + x_arr(randi(numSamples,1));
    z_arr(i) = z;
    y_arr(i) = y;
end

R = var(x_arr);
%Kalman filter
xhat = 0;
Sigma = 1;
x_est_arr = zeros(1,tspan);
for i = 1:tspan
    xhat_pred = xhat;
    Sigma_pred = Sigma + 1;
    S = Sigma_pred + R;
    xhat = xhat_pred + Sigma_pred*(y_arr(i) - xhat_pred)/S;
    Sigma = Sigma_pred - Sigma_pred^2/S;
    x_est_arr(i) = xhat;
end
%particle filter
particalSize = 2000;
zhat = 2*rand(particalSize,1);
x_est_particle = zeros(1,tspan);
for i = 1:tspan
    zhat_pred = zhat + normrnd(0,1,particalSize,1);
    w = cos(y_arr(i)-zhat_pred).^2.*sin(2*(y_arr(i)-zhat_pred)).^2.*normpdf(y_arr(i)-zhat_pred);
    partition = mnrnd(particalSize,w/sum(w));
    zhat = [];
    for j = 1:length(partition)
        num = partition(j);
        for k = 1:num
            zhat = [zhat,zhat_pred(j)];
        end
    end
    x_est_particle (i) = mean(zhat);
end
figure
plot(1:tspan,z_arr,1:tspan,x_est_arr)
hold on
plot(1:tspan,x_est_particle)
xlabel('Time')
legend('True State','Kalman', 'Particle');
figure
plot(1:tspan,z_arr,1:tspan,x_est_arr)
xlabel('Time')
legend('True State','Kalman')