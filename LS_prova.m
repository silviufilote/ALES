clear;
clc;
close all;
rng("default");

%% parameters
theta = 3;
T = 1400;
x_0 = 1;
n = 2500;

%% Settings
std = 0.1;
phi = zeros(n, T);
y = zeros(n, T);

theta_LS_ite = zeros(1,n);
theta_LS_closed = zeros(1,n);

for x = 1:n
    noise_y = 2 * randn(1, T);
    noise_x = std * randn(1, T);
    for t = 1:T
        phi(x,t) = x_0 + noise_x(t);
        y(x,t) = theta * phi(x,t) - theta * noise_x(t) + noise_y(t);
    end
    
    % Theta in forma iterativa
    theta_LS_ite(x) =  pinv(sum(phi(x,:) .* phi(x,:))) * (sum(phi(x,:) .* y(x,:)));
    
    % Theta in forma chiusa
    theta_LS_closed(x) = (phi(x,:) * phi(x,:)') \ (phi(x,:) * y(x,:)');
end



