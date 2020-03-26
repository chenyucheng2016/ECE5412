clc,clear,close all
x0 = normrnd(0,1000,[1,1000]);
T = 1000;
omega_var = 2;
x = x0;
x_next = [];
for i = 1:T
    for j = 1:length(x)
        x_next = [x_next,normrnd(sin(x(j)),omega_var)];
    end
    x = x_next;
    hist(x,50)
end








