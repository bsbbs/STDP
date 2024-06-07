rng(2024);
num_samples = 1600;
% Define the mean and covariance matrix
mean = [0, 0];
covariance = [1, 0.1; 0.1, 1];
covariance = [1, 0.8; 0.8, 1];
data = mvnrnd(mean, covariance, num_samples);
data(data < 0) = 0;
% Separate the data into two vectors
x = data(:, 1);
y = data(:, 2);
plot(x, y, '.');



