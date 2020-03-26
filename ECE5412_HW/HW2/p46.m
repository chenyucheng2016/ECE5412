clc,clear,close all
data = table2array(readtable('dataset'));
X = data(:,2:5);
Y = data(:,6);
var_ests = [];
for i = 1:4
    X_cur = X(:,1:i);
    theta = inv(X'*X)*X'*Y;
    var_est = 1/(length(Y)-i)*(Y-X*theta)'*(Y-X*theta);
    var_ests = [var_ests,var_est]; 
end

plot([0,1,2,3],fliplr(var_ests),'*')
xlabel('# of constraints')
ylabel('Estimated var')