clc,clear,close all
T = 0.1;
A = [1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];
C = [1 0 0 0;0 0 1 0];
f = [T^2/2 0;T 0;0 T^2/2;0 T];
r = [1;2];
z0 = [0;0.1;0;-0.1];
z = z0;
sigma_ob = 1;
sigma_process = 0.1;
Q = diag(sigma_process^2*[1,1,1,1]);
R = diag(sigma_ob^2*[1,1]);
P = eye(4);
z_est0 = zeros(4,1);
z_est = z_est0;
z_array = [];
z_est_array = [];
tspan = 10000;
for i = 1:tspan
    omega = normrnd(0,sigma_process,[4,1]);
    niu = normrnd(0,sigma_ob,[2,1]);
    z = A*z + f*r + omega;
    y = C*z + niu;
    z_est = A*z_est + f*r;
    P = A*P*A'+Q;
    K = P*C'*inv(C*P*C' + R);
    z_est = z_est + K*(y-C*z_est);
    P = (eye(4) - K*C)*P;
    z_array = [z_array,z];
    z_est_array = [z_est_array,z_est];
end
figure(1)
plot(T*(1:tspan),z_array(1,:),T*(1:tspan),z_est_array(1,:))
legend('z(1)','z_{est}(1)')
figure(2)
plot(T*(1:tspan),z_array(2,:),T*(1:tspan),z_est_array(2,:))
legend('z(2)','z_{est}(2)')
figure(3)
plot(T*(1:tspan),z_array(3,:),T*(1:tspan),z_est_array(3,:))
legend('z(3)','z_{est}(3)')
figure(4)
plot(T*(1:tspan),z_array(4,:),T*(1:tspan),z_est_array(4,:))
legend('z(4)','z_{est}(4)')
mean(abs(z_array(1,:)-z_est_array(1,:)))
mean(abs(z_array(2,:)-z_est_array(2,:)))
mean(abs(z_array(3,:)-z_est_array(3,:)))
mean(abs(z_array(4,:)-z_est_array(4,:)))