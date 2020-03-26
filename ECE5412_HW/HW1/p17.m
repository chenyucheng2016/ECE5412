clc,clear,close all
samples = [];
sampleSize = 50000;
for i = 1:sampleSize
    u = rand();
    if u <= (3 - exp(-2))/3
        p = find(mnrnd(1,[1/3,2/3]));
        if p == 1
            x = -log(1 - u)/2;
            if x <= 1
                samples = [samples,x];
            end
        elseif p == 2
            x = u;
            samples = [samples,x];
        end
        
    else
        x = -log(3 - 3*u)/2;
        samples = [samples,x];
    end
end
hist(samples,200);

figure
x_data = 0:0.01:5;
F = [];
for x = x_data
    if x <= 1
        F = [F,(2*exp(-2*x) + 2)/3];
    else
        F = [F,2*exp(-2*x)/3];
    end
end
plot(x_data,F)
