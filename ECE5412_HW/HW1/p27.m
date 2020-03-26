clc,clear,close all
pi0 = [0.1,0.9,0];
A = [0.3,0.3,0.4;0.1,0.9,0;0.1,0.1,0.8];
%% section a
cur_pi = pi0;
samples = [];
for i = 1:500
    samples = [samples,sampleFunc(cur_pi)];
    cur_pi = cur_pi*A;
end
function samp = sampleFunc(cur_pi)
u = rand();
if u < cur_pi(1)
    samp = 1;
elseif u < (cur_pi(1) + cur_pi(2))
    samp = 2;
else
    samp = 3;
end

end
%% section b
proposalDist = [1/3,1/3,1/3];
sampleRej = [];
cur_pi = pi0;
while length(sampleRej) < 500
    c = max(cur_pi)*3;
    x = find(mnrnd(1,proposalDist));
    u = rand();
    if u < cur_pi(x)/(c/3)
        sampleRej = [sampleRej,x];
    end
    cur_pi = cur_pi*A;
end

