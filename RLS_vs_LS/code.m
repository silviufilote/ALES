clear;
clc;
close all;
rng("default");

%% Parameters
T = 500;
b = 3;
u = 1;
theta = b;

%% Initialize vectors 
y = zeros(1, T);
phi = ones(1, T);
noise = randn(1, T);

for t= 1:T
    y(t) = phi(t) * theta + noise(t);
end

%% LS FORM - batch updates
theta_LS = zeros(1,T);
for t = 1:T
    theta_LS(t) = 1/(sum(phi(1:t) .* phi(1:t))) * (sum(phi(1:t) .* y(1:t)));
end

%% Recursive Least Squares (RLS) - form 1

theta_RLS_1 = zeros(1, T);
S = 1;
theta_1 = 0;
for t = 2:T
    S = S + phi(t) * phi(t)';
    epsilon = y(t) - phi(t)' * theta_1;
    K = 1/S * phi(t);
    theta_1 = theta_1 + K * epsilon;
    theta_RLS_1(t) = theta_1;
end

%% Recursive Least Squares (RLS) - form 3
% recursive we update the estimates each time a new data arrives

P = ones(1,T);
theta_RLS_3 = zeros(1, T);
theta_3 = 0;
theta_RLS_3(1,1) = theta_3;
for t = 2:T
    beta = 1 + phi(t)' * P(t - 1) * phi(t);
    P(t) = P(t - 1) - 1/beta * P(t - 1) * phi(t) * phi(t)' * P(t - 1);
    epsilon = y(t) - phi(t)' * theta_3;
    K = P(t) * phi(t);
    theta_3 = theta_3 + K * epsilon;
    theta_RLS_3(t) = theta_3;
end

%% Plotting

figure
hold on 
legend

plot(y,'DisplayName','Noisy data');
p1 = plot(theta_LS,'g','DisplayName','LS');
p2 = plot(theta_RLS_3,'r','DisplayName','RLS3');
% p3 = plot(3 .* ones(1,T), "--", "color", "black", 'HandleVisibility','off');
p3 = yline(3, "--", "color", "black", 'HandleVisibility','off');
p4 = plot(theta_RLS_1,'y','DisplayName','RLS1');
p1.LineWidth = 1;
p2.LineWidth = 2;
p3.LineWidth = 1.5;
p4.LineWidth = 1;
set(legend('Interpreter','Latex'))
xlabel('Time stamp') 
ylabel('$\hat{\theta}$', 'Interpreter','Latex') 
ylim([0, 6]);
xlim([0, 500]);
xticks(0:100:500)