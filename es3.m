clear;
clc;
close all;
rng("default");

%% parameters
theta = 3;
T = 2500;
n = 1200;

%% Sigma_x = 0
x_0 = randn(n, T);
y = zeros(n, T);
theta_LS_1 = zeros(1,n);

for x = 1:n
    noise_y = 2 * randn(1, T);
    for t = 1:T
        y(x,t) = theta * x_0(x,t) + noise_y(t);
    end
    
    % Theta in forma iterativa
    theta_LS_1(x) =  pinv(sum(x_0(x,:) .* x_0(x,:))) * (sum(x_0(x,:) .* y(x,:)));
end


%% Sigma = 0.1
std = 0.1;
x_0 = randn(n, T);
y = zeros(n, T);
theta_LS_2 = zeros(1,n);

for x = 1:n
    noise_y = 2 * randn(1, T);
    noise_x = std * randn(1, T);
    for t = 1:T
        phi(x,t) = x_0(x,t) + noise_x(t);
        % y(x,t) = theta * phi(x,t) - theta * noise_x(t) + noise_y(t);
        y(x,t) = theta * x_0(x,t) + noise_y(t);
    end
    
    % Theta in forma iterativa
    theta_LS_2(x) =  pinv(sum(phi(x,:) .* phi(x,:))) * (sum(phi(x,:) .* y(x,:)));
end
    

%% Sigma = 0.3
std = 0.3;
x_0 = randn(n, T);
y = zeros(n, T);
theta_LS_3 = zeros(1,n);

for x = 1:n
    noise_y = 2 * randn(1, T);
    noise_x = std * randn(1, T);
    for t = 1:T
        phi(x,t) = x_0(x,t) + noise_x(t);
        y(x,t) = theta * phi(x,t) - theta * noise_x(t) + noise_y(t);
        % y(x,t) = theta * x_0(x,t) + noise_y(t);
    end
    
    % Theta in forma iterativa
    theta_LS_3(x) = (1/sum(phi(x,:) .* phi(x,:))) * (sum(phi(x,:) .* y(x,:)));
end

%% Plotting

% First line
tiledlayout(3,1)
nexttile
hold on
title('LS estimate of $\theta$ with $\sigma^2_x = 0$', 'Interpreter', 'latex')
histogram(theta_LS_1, 'DisplayName', 'Histogram of $\hat{\theta}$')
xline(3,'--k', 'LineWidth', 2, 'DisplayName', 'True value of $\theta$');
set(legend('Interpreter','Latex'))
l1 = legend
l1.Location = 'southeast';

ylim([0, 200]);
xticks(0:100:200);
xlim([2, 4]);
xticks(2:0.2:4);


% second line
nexttile
hold on
title('LS estimate of $\theta$ with $\sigma^2_x = 0.1$', 'Interpreter', 'latex')
histogram(theta_LS_2, 'DisplayName', 'Histogram of $\hat{\theta}$')
xline(3,'--k', 'LineWidth', 2, 'DisplayName', 'True value of $\theta$');
set(legend('Interpreter','Latex'))
l2 = legend
l2.Location = 'southeast';

ylim([0, 200]);
xticks(0:100:200);
xlim([2, 4]);
xticks(2:0.2:4);

% third line
nexttile
hold on
title('LS estimate of $\theta$ with $\sigma^2_x = 0.3$', 'Interpreter', 'latex')
histogram(theta_LS_3, 'DisplayName', 'Histogram of $\hat{\theta}$')
xline(3,'--k', 'LineWidth', 2, 'DisplayName', 'True value of $\theta$');
set(legend('Interpreter','Latex'))
l2 = legend
l2.Location = 'southeast';

ylim([0, 200]);
xticks(0:100:200);
xlim([2, 4]);
xticks(2:0.2:4);


