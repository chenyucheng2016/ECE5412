clc,clear,close all
tspan = 1000;
S = [0, 10, 20];
P = [0.1,0.1,0.8;0.3,0.3,0.4;0.2,0.2,0.6];
b0 = [0.3,0.3,0.4];
sigma_sq = 5;
pi_s = b0;
b = b0;
p_error = 0;
ys = [];
for i =1:tspan
    b = b*P;
    x = 10*(find(mnrnd(1,b))-1);
    y = normrnd(x,sqrt(sigma_sq));
    ys = [ys,y];
    pi_s_unnorm = [normpdf(y,0,sqrt(sigma_sq))*pi_s*P(:,1),...
            normpdf(y,10,sqrt(sigma_sq))*pi_s*P(:,2),...
            normpdf(y,20,sqrt(sigma_sq))*pi_s*P(:,3)];
    pi_s =  pi_s_unnorm/sum(pi_s_unnorm);
    p_error = 1 - max(pi_s) + p_error;
end
p_error_tot = 0;
for i = 1:tspan
    %forward
    pi_s = b0;
    for j = 1:i
        pi_s_unnorm = [normpdf(ys(j),0,sqrt(sigma_sq))*pi_s*P(:,1),...
            normpdf(ys(j),10,sqrt(sigma_sq))*pi_s*P(:,2),...
            normpdf(ys(j),20,sqrt(sigma_sq))*pi_s*P(:,3)];
        pi_s =  pi_s_unnorm/sum(pi_s_unnorm);
    end
    %backward
    beta = [1,1,1]';
    for k = tspan:-1:i
        beta = P*diag([normpdf(ys(k),0,sqrt(sigma_sq)),normpdf(ys(k),10,sqrt(sigma_sq)),normpdf(ys(k),20,sqrt(sigma_sq))])*beta;
        beta = beta/sum(beta);
    end
    pi_tot = [pi_s(1)*beta(1),pi_s(2)*beta(2),pi_s(3)*beta(3)];
    pi_tot = pi_tot/sum(pi_tot);
    p_error_tot = 1 - max(pi_tot) + p_error_tot;
end
ave_p_error = p_error/tspan
ave_p_error_tot = p_error_tot/tspan


