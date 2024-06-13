function values = ParetoSequence(num_samples, sigma)
% Define the parameters for the Pareto distribution
alpha = 3; % Shape parameter for the Pareto distribution
scale = 1; % Scale parameter for the Pareto distribution

% Generate a multivariate normal distribution
mean = [0, 0];
covariance = [1, sigma; sigma, 1];
data = mvnrnd(mean, covariance, num_samples);

% Transform the normal distribution to a Pareto distribution
u = normcdf(data); % Convert to uniform distribution
x = scale ./ (1 - u(:, 1)).^(1 / alpha); % Transform to Pareto distribution
y = scale ./ (1 - u(:, 2)).^(1 / alpha); % Transform to Pareto distribution
x = x - 1;
y = y - 1;
values = [x, y];
% % Save the vectors to a .mat file
% save('multivariate_pareto.mat', 'x', 'y');

% % Display the first few rows of the vectors
% disp('First few rows of the multivariate Pareto distributions:');
% disp(table(x(1:5), y(1:5)));

% Plot the data
figure;
scatter(x, y, 'filled');
title('Multivariate Pareto Distributions');
xlabel('X');
ylabel('Y');
grid on;