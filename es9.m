clear;
clc;
rng("default");
close all;

load("data_mimo_kung.mat");

set(0, 'defaultAxesFontSize', 9)
set(0, 'DefaultLineLineWidth', 2);
set(0, 'defaultAxesFontSize', 9)
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultlegendInterpreter','latex')

%% Setup
N = 30;
noise_std = 0.02;
p = 2;                  % n outputs
mu = 3;                 % n inputs
q = 10;                 
d = 20;
time = [0:1:N-1];

A = [0.7 0; 0 0.3];
B = [-0.1 0.2 0.6; 0.9 -0.5 -0.4];
C = [0.5 0.2; 0.6 -0.8];
D = zeros(p, mu);

%% Generate data
e1 = noise_std * randn(N, 2); 
e2 = noise_std * randn(N, 2);
e3 = noise_std * randn(N, 2);

% Deterministic
G0z = ss(A, B, C, D, 1);
impulse_response = impulse(G0z, N - 1);

y_det1 = impulse_response(:, :, 1);
y_det2 = impulse_response(:, :, 2);
y_det3 = impulse_response(:, :, 3);

% Simulated = y(t) + noise
y_sim1 = y1;
y_sim2 = y2;
y_sim3 = y3;

%% henkel matrix

H = zeros(p*(q), mu*(d));
for q_ind = 1:q
    for d_ind = 1:d
        t =  (q_ind + d_ind);
        w = [y_sim1(t,1) y_sim2(t,1) y_sim3(t,1); y_sim1(t,2) y_sim2(t,2) y_sim3(t,2)];
        H(2*(q_ind-1)+1:2*(q_ind-1)+1+1, 3*(d_ind-1)+1:3*(d_ind-1)+1+2) = w;
    end    
end
Hqd = H;

% Hankel matrix
[U,S,V] = svd(Hqd);                                      % SVD on Hankel matrix

n = 2;
Un = U(:, 1:n);
Sn = S(1:n, 1:n);
Vn = V(:, 1:n);

%% Estimates

Oq_hat = Un * sqrt(Sn);
Rd_hat = sqrt(Sn) * Vn';

C_hat = Oq_hat(1:p, :);
B_hat = Rd_hat(:, 1:mu);
D_hat = [y_sim1(1,1) y_sim2(1,1) y_sim3(1,1); y_sim1(1,2) y_sim2(1,2) y_sim3(1,2)];
A_hat = pinv(Oq_hat(1:(q - p), :)) * Oq_hat((p + 1): q, :);

% ss: creates the discrete-time state-space model object 
% tf: creates a continuous-time transfer function model
G_hat = ss(A_hat, B_hat, C_hat, D_hat, 1);
y_sim_det = impulse(G_hat, N - 1);

y_sim_det1 = y_sim_det(:,:,1);
y_sim_det2 = y_sim_det(:,:,2);
y_sim_det3 = y_sim_det(:,:,3);


%% Plotting
figure('Name','Noisy data')
tiledlayout(3,3)
nexttile 
plot(time, y_det1(:,1), 'k--', 'DisplayName', "True impulse response");
hold on;
plot(time, y_sim_det1(:,1), 'b',  'DisplayName', "Estimated impulse response");
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([-0.1,0.2]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 1 to output 1");
set(legend('Interpreter','Latex'))

nexttile 
plot(time, y_det2(:,1), 'k--');
hold on;
plot(time, y_sim_det2(:,1), 'b');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.07]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 2 to output 1");


nexttile 
plot(time, y_det3(:,1), 'k--');
hold on;
plot(time, y_sim_det3(:,1), 'b');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.3]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 3 to output 1");

nexttile 
plot(time, y_det1(:,2), 'k--');
hold on;
plot(time, y_sim_det1(:,2), 'b');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([-1, 0]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 1 to output 2");

nexttile 
plot(time, y_det2(:,2), 'k--');
hold on;
plot(time, y_sim_det2(:,2), 'b');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.6]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 2 to output 2");


nexttile 
plot(time, y_det3(:,2), 'k--');
hold on;
plot(time, y_sim_det3(:,2), 'b');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.8]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 3 -> to output 2");

nexttile 
plot(svd(Hqd), 'b--o');
hold on;
xlabel('Singular value number');
ylabel('singular value');
grid on;
xlim([1,20]);
title("Singular values");

%% Deterministic

y_sim1 = y_det1;
y_sim2 = y_det2;
y_sim3 = y_det3;

H = zeros(p*(q), mu*(d));
for q_ind = 1:q
    for d_ind = 1:d
        t =  (q_ind + d_ind);
        w = [y_sim1(t,1) y_sim2(t,1) y_sim3(t,1); y_sim1(t,2) y_sim2(t,2) y_sim3(t,2)];
        H(2*(q_ind-1)+1:2*(q_ind-1)+1+1, 3*(d_ind-1)+1:3*(d_ind-1)+1+2) = w;
    end    
end
Hqd = H;

% Hankel matrix
[U,S,V] = svd(Hqd);                                      % SVD on Hankel matrix

n = 2;
Un = U(:, 1:n);
Sn = S(1:n, 1:n);
Vn = V(:, 1:n);

Oq_hat = Un * sqrt(Sn);
Rd_hat = sqrt(Sn) * Vn';

C_hat = Oq_hat(1:p, :);
B_hat = Rd_hat(:, 1:mu);
D_hat = [y_sim1(1,1) y_sim2(1,1) y_sim3(1,1); y_sim1(1,2) y_sim2(1,2) y_sim3(1,2)];
A_hat = pinv(Oq_hat(1:(q - p), :)) * Oq_hat((p + 1): q, :);

% ss: creates the discrete-time state-space model object 
% tf: creates a continuous-time transfer function model
G_hat = ss(A_hat, B_hat, C_hat, D_hat, 1);
y_sim_det = impulse(G_hat, N - 1);

y_sim_det1 = y_sim_det(:,:,1);
y_sim_det2 = y_sim_det(:,:,2);
y_sim_det3 = y_sim_det(:,:,3);


%% Plotting
figure('Name','Noiseless data')
tiledlayout(3,3)
nexttile 
plot(time, y_sim_det1(:,1), 'b',  'DisplayName', "Estimated impulse response");
hold on;
plot(time, y_det1(:,1), 'k--', 'DisplayName', "True impulse response");
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([-0.1,0.2]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 1 to output 1");
set(legend('Interpreter','Latex'))

nexttile 
plot(time, y_sim_det2(:,1), 'b');
hold on;
plot(time, y_det2(:,1), 'k--');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.04]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 2 to output 1");


nexttile 
plot(time, y_sim_det3(:,1), 'b');
hold on;
plot(time, y_det3(:,1), 'k--');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.3]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 3 to output 1");

nexttile 
plot(time, y_sim_det1(:,2), 'b');
hold on;
plot(time, y_det1(:,2), 'k--');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([-1, 0]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 1 to output 2");

nexttile 
plot(time, y_sim_det2(:,2), 'b');
hold on;
plot(time, y_det2(:,2), 'k--');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.6]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 2 to output 2");


nexttile 
plot(time, y_sim_det3(:,2), 'b');
hold on;
plot(time, y_det3(:,2), 'k--');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.8]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 3 -> to output 2");

nexttile 
plot(svd(Hqd), 'b--o');
hold on;
xlabel('Singular value number');
ylabel('singular value');
grid on;
xlim([1,20]);
title("Singular values");
