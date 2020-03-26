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

dims = 1:20;
aics = [];
for d = dims
    aic = computeAIC(ys,d,N);
    aics = [aics,aic];
end
figure(1)
plot(aics,'*')
xlabel('Model dimension')
ylabel('AIC score')
bics = [];
for d = dims
    bic = computeBIC(ys,d,N);
    bics = [bics,bic];
end
figure
plot(bics,'*')
xlabel('Model dimension')
ylabel('BIC score')
function aic = computeAIC(ys,d,N)
theta0 = rand(1,d)';
rho = 0.999;
c = 2;
theta = theta0;
P = eye(d);
for i = (d+1):N
    phi = [];
    for j = 1:d
        phi = [phi;ys(i-j)];
    end
    theta = theta + P*phi/(rho/c + phi'*P*phi)*(ys(i) - phi'*theta);
    P = 1/rho*(P - P*(phi*phi')*P/(rho/c + phi'*P*phi));
end
ysim = ys(1:d);
for i = (d+1):N
    phi = [];
    for j = 1:d
        phi = [phi;ys(i-j)];
    end
    y = theta'*phi;
    ysim = [ysim,y];
end
aic = norm(ys-ysim) + 2*d;
end


function bic = computeBIC(ys,d,N)
theta0 = rand(1,d)';
rho = 0.999;
c = 2;
theta = theta0;
P = eye(d);
for i = (d+1):N
    phi = [];
    for j = 1:d
        phi = [phi;ys(i-j)];
    end
    theta = theta + P*phi/(rho/c + phi'*P*phi)*(ys(i) - phi'*theta);
    P = 1/rho*(P - P*(phi*phi')*P/(rho/c + phi'*P*phi));
end
ysim = ys(1:d);
for i = (d+1):N
    phi = [];
    for j = 1:d
        phi = [phi;ys(i-j)];
    end
    y = theta'*phi;
    ysim = [ysim,y];
end
bic = log(N)*d + 2*(N/2*log(2*pi)+1/2*norm(ys-ysim)); 
end