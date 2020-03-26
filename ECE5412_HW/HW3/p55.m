%problem 55
clear,clc,close all
data = importdata('australian');
data = data.textdata;
[numData, numFeature] = size(data);
dataset = zeros(numData,numFeature);
for i = 1:numData
    for j = 1:numFeature
        cell_data = str2num(data{i,j});
        dataset(i,j) = cell_data(1);
    end 
end
randIndices = randperm(numData);
trainingNum = round(numData*0.9);
trainingData  = dataset(randIndices(1:trainingNum),:);
testData = dataset(randIndices(trainingNum+1:end),:);
MdLinear = fitcdiscr(trainingData(:,2:end),trainingData(:,1));
predictedClass_LDA = predict(MdLinear, testData(:,2:end));
predictedClass_LDA = max(predictedClass_LDA,0);
x0 = rand(1,13)/100;
trainingFeature = trainingData(:,2:end);
trainingLabel = max(trainingData(:,1),0);
options = optimset('PlotFcns',@optimplotfval,'MaxIter',100000,'MaxFunEvals',100000);
[x,val,etflag] = fminunc(@(x)costFcn(x,trainingFeature,trainingLabel),x0,options);
predicted_LR = 1./(1+exp(-testData(:,2:end)*x'));
for i = 1:length(predicted_LR)
    if predicted_LR(i) >= 0.5
        predicted_LR(i) = 1;
    else
        predicted_LR(i) = 0;
    end
end
testLabel = max(testData(:,1),0);
LDA_rate= 0;
LR_rate = 0;
for i = 1:length(testLabel)
    if abs(predicted_LR(i)-testLabel(i)) < 0.01
        LR_rate = LR_rate + 1;
    end
    if abs(predictedClass_LDA(i) - testLabel(i)) < 0.01
        LDA_rate = LDA_rate + 1;
    end
end


function cost = costFcn(x,trainingData,trainingLabel)
cost = 0;
numData = size(trainingData,1);
for i = 1:numData
    dataline = trainingData(i,:);
    sumLine = sum(dataline.*x);
    y = trainingLabel(i);
    p = 1/(1+exp(-sumLine));
    cost = cost - (y*log(p) + (1-y)*log(1-p));
end
end



