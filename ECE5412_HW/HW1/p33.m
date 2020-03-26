clc,clear,close all
N = 1000;
Ms = [100,500,1000,2000,5000,8000,10000] ;
L = 100;
A = rand(N);
pik = rand([N,1]);
pik = pik/sum(pik);
%Use rejection sampler
proposalDist = ones(1,N)/N;
c = max(pik)/(proposalDist(1));
At = A';
piKp1_est = zeros(N,1);
tic 
piKp1 = At*pik;
t1 = toc;
mses = [];
vars = [];
times = [];
for M = Ms
    errs = [];
    piKp1_estm = [];
    ts = [];
    for i = 1:L
        num = 0;
        tic
        while (num < M)
            x = find(mnrnd(1,proposalDist));
            u = rand();
            if u < pik(x)/(c/N)
                %accept
                piKp1_est = piKp1_est + At(:,x);
                num = num + 1;
            end
        end
        piKp1_est = piKp1_est/M;
        t = toc;
        ts = [ts,t];
        piKp1_estm = [piKp1_estm,piKp1_est];
        err = norm(piKp1 - piKp1_est);
        errs = [errs,err];
    end
    times = [times,mean(ts)];
    meanEst = mean(piKp1_estm,2);
    var = trace((piKp1_estm - meanEst)'*(piKp1_estm - meanEst))/L;
    mse = sum(errs)/L;
    mses = [mses,mse];
    vars = [vars,var];
end
figure
plot(Ms,vars,'*');
xlabel('# of iterations')
ylabel('estimator variance')
figure
plt1 = plot(Ms,times,'b*');
hold on
plt2 = line([0,10000],[t1,t1],'LineWidth',1.5','Color','red');
legend([plt1,plt2],{'Randomized method','Direct multiplication'})
xlabel('# of iterations')
ylabel('Running time [sec]')
