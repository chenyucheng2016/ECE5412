clc,clear,close all
%question 2a
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

step_cdf = 0.0001;
lower_bound = -4;
upper_bound = 4;
cdf = [];
cum_den = 0;
for i = lower_bound:step_cdf:upper_bound
    cum_den = cum_den + cos(i)^2*sin(2*i)^2*normpdf(i)*step_cdf;
    cdf = [cdf,cum_den];
end
cdf = cdf/max(cdf);
figure
[f,x_h] = ecdf(x_arr);
plot(x_h,f,'b','LineWidth',2);
hold on
plot(lower_bound:step_cdf:upper_bound,cdf,'r--','LineWidth',2)
legend('Empirical CDF','True CDF');