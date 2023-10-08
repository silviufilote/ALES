clear;
clc;
close all;
rng("default");

%% settings signals
T = 1000;
u_1 = 10 * randn(1, T);
u_2 = 10 * randn(1, T); 
e = randn(1, T); 
theta_2 = ones(1,T);

theta_1 = zeros(1,T);
for t =1:T
    if t < 400
        theta_1(t) = 1;
    else
        theta_1(t) = sin(0.01 * pi * t + pi/2);
    end
end

y = zeros(1,T);
for t = 1:T
    y(t) = theta_1(t) * u_1(t) + theta_2(t) * u_2(t) + e(t);
end


%% Estimate parameters using RLS form 3 - with and without forgetting factor 

phi = zeros(2,T);
for t = 1:T
    phi(:,t) = [u_1(t), u_2(t)]';
end

theta = zeros(2,T);
for t = 1:T
    theta(:,t) = [theta_1(t), theta_2(t)]';
end

theta_RLS = zeros(2, T);
P = diag([10^6 10^6]);
theta_0 = [0, 0]';
for t = 1:T
    beta = 1 + phi(:,t)'  * P * phi(:,t);
    P = P - 1/beta * P * phi(:,t) * phi(:,t)' * P;
    epsilon = y(t) - phi(:,t)' * theta_0;
    K = P * phi(:,t);
    theta_0 = theta_0 + K * epsilon;
    theta_RLS(:,t) = theta_0;
end

theta_RLS_2 = zeros(2, T);
P = diag([10^6 10^6]);
theta_0 = [0, 0]';
mu = 0.95;
for t = 1:T
    beta = mu + phi(:,t)'  * P * phi(:,t);
    P = 1/mu * (P - 1/beta * P * phi(:,t) * phi(:,t)' * P);
    epsilon = y(t) - phi(:,t)' * theta_0;
    K = P * phi(:,t);
    theta_0 = theta_0 + K * epsilon;
    theta_RLS_2(:,t) = theta_0;
end

theta_RLS_3 = zeros(2, T);
P = diag([10^6 10^6]);
theta_0 = [0, 0]';
mu = 0.75;
for t = 1:T
    beta = mu + phi(:,t)'  * P * phi(:,t);
    P = 1/mu * (P - 1/beta * P * phi(:,t) * phi(:,t)' * P);
    epsilon = y(t) - phi(:,t)' * theta_0;
    K = P * phi(:,t);
    theta_0 = theta_0 + K * epsilon;
    theta_RLS_3(:,t) = theta_0;
end

%% LMS 
theta_LMS = zeros(2, T);
theta_0 = [0, 0]';
gamma = 0.002;
for t = 1:T
    epsilon = y(t) - phi(:,t)' * theta_0;
    K = gamma * phi(:,t);
    theta_0 = theta_0 + K * epsilon;
    theta_LMS(:,t) = theta_0;
end

%% Plotting 

% First line
tiledlayout(4,2)
nexttile
hold on
plot(theta_1, "--","color", "black", 'DisplayName', 'True value $\hat{\theta}_1$');
plot(theta_RLS(1, :), 'blue', 'DisplayName', 'Estimated value $\hat{\theta}_1$');
xlabel('Time stamp') 
ylabel('RLS with \mu = 1') 
set(legend('Interpreter','Latex'))
l1 = legend
l1.Location = 'southwest';

ylim([-1, 1.5]);
xlim([0, 1000]);
xticks(0:200:1000)

nexttile
hold on
plot(theta_2, "--","color", "black", 'DisplayName', 'True value $\hat{\theta}_2$');
plot(theta_RLS(2, :), 'blue', 'DisplayName', 'Estimated value $\hat{\theta}_2$');
xlabel('Time stamp') 
ylabel('RLS with \mu = 1') 
set(legend('Interpreter','Latex'))
l2 = legend
l2.Location = 'southeast';

ylim([-1, 1.5]);
xlim([0, 1000]);
xticks(0:200:1000)

% 2nd row

nexttile
hold on
plot(theta_1, "--","color", "black", 'DisplayName', 'True value $\hat{\theta}_1$');
plot(theta_RLS_2(1, :), 'blue', 'DisplayName', 'Estimated value $\hat{\theta}_1$');
xlabel('Time stamp') 
ylabel('RLS with \mu = 0.95') 

ylim([-1, 1.5]);
xlim([0, 1000]);
xticks(0:200:1000)

nexttile
hold on
plot(theta_2, "--","color", "black", 'DisplayName', 'True value $\hat{\theta}_2$');
plot(theta_RLS_2(2, :), 'blue', 'DisplayName', 'Estimated value $\hat{\theta}_2$');
xlabel('Time stamp') 
ylabel('RLS with \mu = 0.95') 

ylim([-1, 1.5]);
xlim([0, 1000]);
xticks(0:200:1000)

% 3rd row 

nexttile
hold on
plot(theta_1, "--","color", "black", 'DisplayName', 'True value $\hat{\theta}_1$');
plot(theta_RLS_3(1, :), 'blue', 'DisplayName', 'Estimated value $\hat{\theta}_1$');
xlabel('Time stamp') 
ylabel('RLS with \mu = 0.75') 

ylim([-1, 1.5]);
xlim([0, 1000]);
xticks(0:200:1000)

nexttile
hold on
plot(theta_2, "--","color", "black", 'DisplayName', 'True value $\hat{\theta}_2$');
plot(theta_RLS_3(2, :), 'blue', 'DisplayName', 'Estimated value $\hat{\theta}_2$');
xlabel('Time stamp') 
ylabel('RLS with \mu = 0.75') 

ylim([-1, 1.5]);
xlim([0, 1000]);
xticks(0:200:1000)

% 4th row

nexttile
hold on
plot(theta_1, "--","color", "black", 'DisplayName', 'True value $\hat{\theta}_2$');
plot(theta_LMS(1, :), 'blue', 'DisplayName', 'Estimated value $\hat{\theta}_2$');
xlabel('Time stamp') 
ylabel('LMS with \gamma = 0.002') 

ylim([-1, 1.5]);
xlim([0, 1000]);
xticks(0:200:1000)


nexttile
hold on
plot(theta_2, "--","color", "black", 'DisplayName', 'True value $\hat{\theta}_2$');
plot(theta_LMS(2, :), 'blue', 'DisplayName', 'Estimated value $\hat{\theta}_2$');
xlabel('Time stamp') 
ylabel('LMS with \gamma = 0.002') 

ylim([-1, 1.5]);
xlim([0, 1000]);
xticks(0:200:1000)