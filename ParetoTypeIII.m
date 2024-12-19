function samples = pareto_type_III(N, xm, alpha, beta)
    % Simulate N samples from the Pareto Type III distribution
    % Parameters:
    % N     : Number of samples
    % xm    : Scale parameter
    % alpha : Shape parameter
    % beta  : Location parameter
    
    % Step 1: Generate uniform random numbers
    U = rand(N, 1);
    
    % Step 2: Apply inverse CDF transformation
    samples = beta + xm * ((1 - U).^(-1/alpha) - 1);
end

% Example usage
N = 10000; % Number of samples
xm = 2;    % Scale parameter
alpha = 3; % Shape parameter
beta = 0;  % Location parameter

% Generate Pareto Type III samples
data = pareto_type_III(N, xm, alpha, beta);

% Plot histogram of the simulated data
figure;
histogram(data, 50, 'Normalization', 'pdf');
xlabel('Value');
ylabel('Density');
title('Pareto Type III Distribution Simulation');