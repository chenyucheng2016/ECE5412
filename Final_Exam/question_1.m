%Question 1
clc,clear,close all
P = [0.8,0.2;0.2,0.8];
%P = [0.1, 0.9;0.3,0.7];
pi0 = [0.5,0.5];
tspan = 1000;
num_sim = 50;
A = 0.7;
omega = 5;
errors = zeros(1,num_sim);
errors_filter = zeros(1,num_sim);
for j = 1:num_sim
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
        y(i) = s + randn(1) + A*sin(omega*i);
        pi_d = pi_d*P;
        pi_d = pi_d / sum(pi_d);
    end
    
    p_nxt = 0.5;
    s_est = zeros(1,tspan);
    s_est_filter = zeros(1,tspan);
    backward_prob = zeros(1,tspan);
    forward_prob = zeros(1,tspan);
    gamma_prob = zeros(1,tspan);
    backward_prob(end) = 0.5;
    p_back = [1;1];
    %backward algorithm
    for k = tspan:-1:1
        %B = diag([normpdf(y(k+1),-1 + A*sin(omega*(k+1)),1),normpdf(y(k+1),1 + A*sin(omega*(k+1)),1)]);
        B = diag([normpdf(y(k),-1 + A*sin(omega*(k)),1),normpdf(y(k),1 + A*sin(omega*(k)),1)]);
        p_back = P*B*p_back;
        p_back = p_back/sum(p_back);
        backward_prob(k) = p_back(1);
    end
    %implementation of filter
    for i = 1:tspan
        %forward filter
        p_nxt = normpdf(y(i),-1 + A*sin(omega*i),1)*(p_nxt * P(1,1) + (1 - p_nxt)* P(2,1))/(normpdf(y(i),-1 + A*sin(omega*i),1)*(p_nxt * P(1,1) + (1 - p_nxt)* P(2,1)) + ...
            normpdf(y(i),1 + A*sin(omega*i),1)*(p_nxt * P(1,2) + (1 - p_nxt)* P(2,2)));
        forward_prob(i) = p_nxt;
        s_est_filter(i) = -p_nxt + (1 - p_nxt);
        gamma1 = p_nxt * backward_prob(i);
        gamma2 = (1-p_nxt) * (1 - backward_prob(i));
        gammas = [gamma1, gamma2]/(gamma1 + gamma2);
        gamma_prob(i) = gammas(1);
        s_est(i) = -gammas(1) + gammas(2);
    end
    errors(j) = norm(s_true - s_est)^2/tspan;
    errors_filter(j) = norm(s_true - s_est_filter)^2/tspan;
end
probs = [forward_prob;backward_prob;gamma_prob;s_true];
figure
plot(1:num_sim,errors,'*')
xlabel('Sim trials')
ylabel('Error Smoother')
figure
plot(1:num_sim,errors_filter,'*')
xlabel('Sim trials')
ylabel('Error filter')
mean(errors)
mean(errors_filter)
% plot(1:tspan,s_true,'+',1:tspan,s_est,'*');
% title('True state vs Estimated State');
% legend('True State','Condition Mean estimate');
% xlabel('Time Steps')
