clc,clear,close all
table = xlsread('AAPL.csv');
price = table(:,5);
logR = [];
for i = 1:size(price,1)-1
    invest = price(i);
    outcome = price(i+1);
    logR = [logR,log(invest/outcome)];
end
h = histogram(logR,50);
p = histcounts(logR,50,'Normalization','pdf');
figure
binCenters = h.BinEdges + (h.BinWidth/2);
plt1 = plot(binCenters(1:end-1), p, 'r-');
sampleMean = mean(logR);
sampleVar = std(logR);
x = -0.2:0.001:0.2;
y = normpdf(x,sampleMean,sampleVar);
ycdf = normcdf(x,sampleMean,sampleVar);
hold on
plt2 = plot(x,y);
legend([plt1,plt2],{'Data Distribution','Normal'})
edges = h.BinEdges;
probMass = [];
pdf = [x;y];
for i = 1:length(edges)-1
    curEdge = edges(i);
    nextEdge = edges(i+1);
    probMass = [probMass,computecdf(pdf,nextEdge)-computecdf(pdf,curEdge)]; 
end
ExpectedCounts = probMass*length(price);
actualCounts = histcounts(logR,50);
chi_square = sum((((ExpectedCounts(1:46) - actualCounts(1:46)).^2)./ExpectedCounts(1:46)));
function comprob = computecdf(pdf,zeta)
xs = pdf(1,:);
comprob = 0;
for i = 1:length(xs)
    x = xs(i);
    if x < zeta
        comprob = comprob + pdf(2,i)*0.0066;
    else
        break;
    end
end
end