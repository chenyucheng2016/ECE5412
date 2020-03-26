clc,clear,close all
a1 = 0.5;
a2 = -0.5;
N = 10000;
y = [];
y0 = 1;
y1 = 4;
ys = [y0,y1];
for i = 3:N
    y = a1*ys(i-1)+a2*ys(i-2) + randn(1);
    ys = [ys,y];
end

theta0 = [-50,20]';
rho = 0.999;
c = 2;
theta = theta0;
P = eye(2);
errs = [];
for i = 3:N
    phi = [ys(i-1);ys(i-2)];
    theta = theta + P*phi/(rho/c + phi'*P*phi)*(ys(i) - phi'*theta);
    P = 1/rho*(P - P*(phi*phi')*P/(rho/c + phi'*P*phi));
    err = norm([theta(1)-a1,theta(2)-a2]);
    errs = [errs,err];
end
plot(errs,'*')
xlabel('Iterations')
ylabel('Error')