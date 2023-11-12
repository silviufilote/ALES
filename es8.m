clear;
clc;
rng("default");
close all;

load("data_siso_kung.mat");

set(0, 'defaultAxesFontSize', 9)
set(0, 'DefaultLineLineWidth', 2);
set(0, 'defaultAxesFontSize', 9)
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultlegendInterpreter','latex')

%% Setup
N = 30;
q = 10;
d = 20;
time = [0:1:N-1];

a1 = 0.2;
a2 = 0.35;
b1 = 1;
b2 = -5;
std = 0.2;

z = tf('z');                                            % creates special variable z that you can use in a rational expression to create a discrete-time transfer function model
G0z = (b1*z^-1 + b2*z^-2)/(1 - a1*z^-1 - a2*z^-2);      % no time sample
H0z = 1/(1 - a1*z^-1 - a2*z^-2);                        % no time sample

%% Generate data
e = std * randn(1, N);

% computes the impulse response y of the dynamic system sys.
y_det = impulse(G0z, N - 1);    % noiseless data
e_sim = lsim(H0z, e);           % simulate noise signal
y_sim = y;                      % data loaded from mat
y_sim2 = y_det + e_sim;

%% Hankel matrix Hqd + SVD                                           
Hqd = hankel(y_sim(2:end,:)); 
Hqd = Hqd(1:q,1:d);

% Hankel matrix
[U,S,V] = svd(Hqd);                                      % SVD on Hankel matrix

n = 2;
Un = U(:, 1:n);
Sn = S(1:n, 1:n);
Vn = V(:, 1:n);

%% Estimates

Oq_hat = Un * sqrt(Sn);
Rd_hat = sqrt(Sn) * Vn';

C_hat = Oq_hat(1, :);
B_hat = Rd_hat(:, 1);
D_hat = y_sim(1);
A_hat = pinv(Oq_hat(1:q - 1, :)) * Oq_hat(2:q,:);


% ss: creates the discrete-time state-space model object 
% tf: creates a continuous-time transfer function model
G_hat = ss(A_hat, B_hat, C_hat, D_hat, 1);
y_sim_hat = impulse(G_hat, N - 1);

%% Plotting
figure('Name','Noisy data')
tiledlayout(1,3)
nexttile 
plot(time, y_sim, 'b', 'DisplayName', "Noisy output");
hold on;
plot(time, y_det, 'k--', 'DisplayName', "Noiseless output");
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([-6,2]); yticks(-6:2:2);
xlim([0,30]); xticks(0:5:time(end));
title("Impulse response");
set(legend('Interpreter','Latex'))

nexttile 
plot(svd(Hqd), 'b--o');
hold on;
xlabel('Singular value number');
ylabel('singular value');
grid on;
xlim([1,10]);
title("Singular values");


nexttile 
plot(time, y_sim_hat, 'b', 'DisplayName', "Estimated impulse response");
hold on;
plot(time, y_det, 'k--', 'DisplayName', "True impulse response");
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([-6,2]); yticks(-6:2:2);
xlim([0,30]); xticks(0:5:time(end));
title("Impulse response");
set(legend('Interpreter','Latex'))

fprintf('True poles: %d, %d \n', pole(minreal(G0z)));
fprintf('Estimated poles: %d, %d \n', pole(minreal(G_hat)));

fprintf('True zero: %d, %d \n', zero(minreal(G0z)));
fprintf('Estimated zero: %d, %d \n', zero(minreal(G_hat)));


%% Deterministic 

y_sim = y_det;

Hqd = hankel(y_sim(2:end,:)); 
Hqd = Hqd(1:q,1:d);

[U,S,V] = svd(Hqd);                                      

n = 2;
Un = U(:, 1:n);
Sn = S(1:n, 1:n);
Vn = V(:, 1:n);

Oq_hat = Un * sqrt(Sn);
Rd_hat = sqrt(Sn) * Vn';

C_hat = Oq_hat(1, :);
B_hat = Rd_hat(:, 1);
D_hat = y_sim(1);
A_hat = pinv(Oq_hat(1:q - 1, :)) * Oq_hat(2:q,:);


% ss: creates the discrete-time state-space model object 
% tf: creates a continuous-time transfer function model
G_hat = ss(A_hat, B_hat, C_hat, D_hat, 1);
y_sim_hat = impulse(G_hat, N - 1);

%% Plotting
figure('Name','Noiseless data')
title("Noise example")
tiledlayout(1,2)

nexttile 
plot(svd(Hqd), 'b--o');
hold on;
xlabel('Singular value number');
ylabel('singular value');
grid on;
xlim([1,10]);
title("Singular values");


nexttile 
plot(time, y_sim_hat, 'b', 'DisplayName', "Estimated impulse response");
hold on;
plot(time, y_det, 'k--', 'DisplayName', "True impulse response");
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([-6,2]); yticks(-6:2:2);
xlim([0,30]); xticks(0:5:time(end));
title("Impulse response");
set(legend('Interpreter','Latex'))


fprintf('True poles: %d, %d \n', pole(minreal(G0z)));
fprintf('Estimated poles: %d, %d \n', pole(minreal(G_hat)));

fprintf('True zero: %d, %d \n', zero(minreal(G0z)));
fprintf('Estimated zero: %d, %d \n', zero(minreal(G_hat)));