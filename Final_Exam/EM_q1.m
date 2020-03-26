%question 1 EM
clc,clear,close
P = [0.8,0.2;0.2,0.8];
pi0 = [0.1,0.9];
tspan = 1000;
num_sim = 50;
A = 0.5;
omega = 5;
pi_d = pi0;
y = zeros(1,tspan);
s_true = zeros(1,tspan);
%genenrate true state sequence and observation vector (y)
for i = 1:tspan
    s = find(mnrnd(1,pi_d)) - 2;
    if s == 0
        s = 1;
    end
    s_true(i) = s;
    y(i) = s + randn + A*sin(omega*i);
    pi_d = pi_d*P;
    pi_d = pi_d / sum(pi_d);
end
%EM algorithm
A = 10;
num_iter = 50;
Logliks = zeros(1,num_iter);
for i = 1:num_iter%iterate 100 times
    % E step
    [gamma_marginal, gamma_joint] = computePosterior(y, A, tspan);
    Logliks(i) = evalLoglik(gamma_marginal, gamma_joint, A, omega, tspan, y);
    % M step
    A = sum((y - (-gamma_marginal + (1 - gamma_marginal))).*sin(omega*(1:tspan)))/sum(sin(omega*(1:tspan)).*sin(omega*(1:tspan)));
end
plot(1:num_iter, Logliks)
ylabel('Log Likelihood')
xlabel('# iter')
function [gamma_marginal, gamma_joint] = computePosterior(y, A, tspan)
forward_prob = zeros(1,tspan);
backward_prob = zeros(1,tspan);
gamma_marginal = zeros(1,tspan);%prob of being -1
gamma_joint = zeros(4,tspan-1);%(-1,-1), (-1,1), (1,-1), (1,1)
P = [0.8,0.2;0.2,0.8];
omega = 5;
p_back = [1;1];
%backward algorithm
for k = tspan:-1:1
    B = diag([normpdf(y(k),-1 + A*sin(omega*k),1),normpdf(y(k),1 + A*sin(omega*k),1)]);
    p_back = P*B*p_back;
    p_back = p_back/sum(p_back);
    backward_prob(k) = p_back(1);
end
p_nxt = 0.5;
for i = 1:tspan
    %forward filter
    p_nxt = normpdf(y(i),-1 + A*sin(omega*i),1)*(p_nxt * P(1,1) + (1 - p_nxt)* P(2,1))/(normpdf(y(i),-1 + A*sin(omega*i),1)*(p_nxt * P(1,1) + (1 - p_nxt)* P(2,1)) + ...
        normpdf(y(i),1 + A*sin(omega*i),1)*(p_nxt * P(1,2) + (1 - p_nxt)* P(2,2)));
    gamma1 = p_nxt*backward_prob(i);
    gamma2 = (1 - p_nxt) * (1 - backward_prob(i));
    gammas = [gamma1, gamma2]/(gamma1 + gamma2);
    gamma_marginal(i) = gammas(1);
    forward_prob(i) = p_nxt;
end

for i = 1:tspan-1
    %s_{i-1} = -1, s_i = -1
    gamma_joint(1,i) = forward_prob(i)*P(1,1)*normpdf(y(i+1),-1 + A*sin(omega*(i+1)),1)*backward_prob(i+1);
    %s_{i-1} = -1, s_i = 1
    gamma_joint(2,i) = forward_prob(i)*P(1,2)*normpdf(y(i+1),1 + A*sin(omega*(i+1)),1)*(1 - backward_prob(i+1));
    %s_{i-1} = 1, s_i = -1
    gamma_joint(3,i) = (1 - forward_prob(i))*P(2,1)*normpdf(y(i+1),-1 + A*sin(omega*(i+1)),1)* backward_prob(i+1);
    %s_{i-1} = 1, s_i = 1
    gamma_joint(4,i) = (1 - forward_prob(i))*P(2,2)*normpdf(y(i+1),1 + A*sin(omega*(i+1)),1)* (1 - backward_prob(i+1));
    gamma_joint(:,i) = gamma_joint(:,i)/sum(gamma_joint(:,i));
end

end

function Loglik = evalLoglik(gamma_marginal, gamma_joint, A, omega, tspan, y)
P = [0.8,0.2;0.2,0.8];
Loglik = -sum(gamma_marginal/2.*(y + 1 - A*sin(omega*(1:tspan))).^2) - sum((1 - gamma_marginal)/2.*(y - 1 - A*sin(omega*(1:tspan))).^2);

for i = 1:size(gamma_joint,1)
    gamma = gamma_joint(:,i);
    Loglik = Loglik + gamma(1)*log(P(1,1)) + gamma(2)*log(P(1,2)) + gamma(3)*log(P(2,1)) + gamma(4)*log(P(2,2));
end
end

