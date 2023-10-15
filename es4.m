clear;
clc;
rng("default");
close all;

%% setup
T = 2500;
n = 1000;
y_1 = zeros(n,T);
y_2 = zeros(n,T);
theta_0 = [0.2, -0.5, 1, 0.5, 0.75]';

theta_LS = ones(3,n);
theta_IV = ones(3,n);


for x = 1:n
    noise_1 = randn(1, T);
    noise_2 = randn(1, T);
    u = randn(1, T);
    
    phi = ones(3,T);
    z = ones(3,T);
    
    for t = 3:T      
         y_1(x,t) = theta_0(1) * y_1(x,t - 1) + theta_0(2) * y_1(x,t - 2) + theta_0(3) *noise_1(t) + theta_0(4) * noise_1(t - 1) + theta_0(5) * noise_1(t - 2) + theta_0(3) * u(t - 1);
         y_2(x,t) = theta_0(1) * y_2(x,t - 1) + theta_0(2) * y_2(x,t - 2) + theta_0(3) *noise_2(t) + theta_0(4) * noise_2(t - 1) + theta_0(5) * noise_2(t - 2) + theta_0(3) * u(t - 1);

         phi(:,t) = [y_1(x,t - 1), y_1(x,t - 2), u(t - 1)]';
         z(:,t) = [y_2(x,t - 1), y_2(x,t - 2), u(t - 1)]';
    end

     % theta_LS(:,x) =  (1\(sum(phi(:,:) * phi(:,:)')) * (sum(phi(:,:) * y_1(x,:)')))';
     % theta_IV(:,x) =  (1\(sum(z(:,:) * z(:,:)')) * (sum(z(:,:) * y_2(x,:)')))';
    
    theta_LS(:,x) = (phi * phi') \ (phi * y_1(x,:)');
    theta_IV(:,x) = (z * phi') \ (z * y_1(x,:)');
end


%% Plotting

% I row
tiledlayout(3,2)
nexttile
hold on
title('LS estimate of $\theta_1$', 'Interpreter', 'latex')
histogram(theta_LS(1,:), 'DisplayName', 'Histogram of $\hat{\theta}_1$')
xline(theta_0(1),'--k', 'LineWidth', 2, 'DisplayName', 'True value of $\theta_1$');
set(legend('Interpreter','Latex'))
l1 = legend
l1.Location = 'southwest';

ylim([0, 200]);
xlim([-0.6,1.1]);
xticks(-0.6:0.2:1.1);


nexttile
hold on
title('IV estimate of $\theta_1$', 'Interpreter', 'latex')
histogram(theta_IV(1,:), 'DisplayName', 'Histogram of $\hat{\theta}_1$')
xline(theta_0(1),'--k', 'LineWidth', 2, 'DisplayName', 'True value of $\theta_1$');
set(legend('Interpreter','Latex'))
l1 = legend
l1.Location = 'southwest';

ylim([0, 200]);
xlim([-0.6,1.1]);
xticks(-0.6:0.2:1.1);


% II row

nexttile
hold on
title('LS estimate of $\theta_2$', 'Interpreter', 'latex')
histogram(theta_LS(2,:), 'DisplayName', 'Histogram of $\hat{\theta}_2$')
xline(theta_0(2),'--k', 'LineWidth', 2, 'DisplayName', 'True value of $\theta_2$');
set(legend('Interpreter','Latex'))
l1 = legend
l1.Location = 'southeast';

ylim([0, 200]);
xlim([-0.6,1.1]);
xticks(-0.6:0.2:1.1);


nexttile
hold on
title('IV estimate of $\theta_2$', 'Interpreter', 'latex')
histogram(theta_IV(2,:), 'DisplayName', 'Histogram of $\hat{\theta}_2$')
xline(theta_0(2),'--k', 'LineWidth', 2, 'DisplayName', 'True value of $\theta_2$');
set(legend('Interpreter','Latex'))
l1 = legend
l1.Location = 'southeast';

ylim([0, 200]);
xlim([-0.6,1.1]);
xticks(-0.6:0.2:1.1);


% III row

nexttile
hold on
title('LS estimate of $\theta_3$', 'Interpreter', 'latex')
histogram(theta_LS(3,:), 'DisplayName', 'Histogram of $\hat{\theta}_3$')
xline(theta_0(3),'--k', 'LineWidth', 2, 'DisplayName', 'True value of $\theta_3$');
set(legend('Interpreter','Latex'))
l1 = legend
l1.Location = 'southwest';

ylim([0, 200]);
xlim([-0.6,1.1]);
xticks(-0.6:0.2:1.1);


nexttile
hold on
title('IV estimate of $\theta_3$', 'Interpreter', 'latex')
histogram(theta_IV(3,:), 'DisplayName', 'Histogram of $\hat{\theta}_3$')
xline(theta_0(3),'--k', 'LineWidth', 2, 'DisplayName', 'True value of $\theta_3$');
set(legend('Interpreter','Latex'))
l1 = legend
l1.Location = 'southwest';

ylim([0, 200]);
xlim([-0.6,1.1]);
xticks(-0.6:0.2:1.1);


