clear all; clc;

input_raw = load('input_data_numeric.csv');
x = input_raw(:,2:end);
w = rand(4,1)*0.4 - 0.2;
theta = rand(1)*2 - 1;
eta = 0.02;

delta = zeros(length(A),1);
A = [1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1];
for mu = 1:length(A)
    V = x(mu,:);
    O(mu) = tanh(0.5*w'*V' - theta);
end
delta(length(A)) = 
for mu = length(A):2
    