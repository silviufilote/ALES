clear;
clc;
rng("default");
close all;

%% Setup
N = 10000;
SNR = 100;
TS = 0.05;
fs = 1/TS; 
pulse = 1e0:(1e-2):(fs*pi); 

% where z is the variable in the Laplace domain
G0z = tf([0 0 0 0.28261 0.50666],[1 -1.41833 1.58939 -1.31608 0.88642],'Ts', TS, 'variable', 'z^-1');
R0z = tf([0.32905 -0.59771 0.70728 -0.64010 0.46499 -0.11769], [1 -1], 'Ts', TS, 'variable', 'z^-1');
H0z = tf([1 0.134], [1 -0.789], 'Ts', TS, 'variable', 'z^-1');
S0z = 1/(1 + G0z*R0z);

%% RPBS
range = [-2 2];
band = [0 0.1];   
step = (0:TS:(N-1)*TS);

r2 = idinput(N, 'prbs', band, range);
r = lsim(R0z, r2, step);

%% Simulate true signal 

Lz = feedback(R0z*G0z, 1);
% Lz = G0z * S0z;
y0_det = lsim(Lz, r2, step);

% lambda noise
var_v = var(y0_det)/SNR;
e = randn(N,1);
H0z_var = lsim(H0z, e, step);

% now v(t) = H0z * e(t) unitary variance noise e(t)
% I divide in this way I can get the actualy noise variance 
var_e = var_v / var(H0z_var);

e0 = sqrt(var_e) * randn(N,1);
y0_sim = lsim(G0z*S0z, r) + lsim(S0z*H0z, e0);
u0 = lsim(S0z, r) - lsim(R0z*S0z*H0z, e0);
data = iddata(y0_sim(500:N), u0(500:N), TS, 'Domain', 'Time');

%% direct identification method -> PEM estimation (open loop)

% SNR = 100
bjM = bj(data,[2 1 1 4 3]);
oeM = oe(data, [2 4 3]);

G_bj1 = tf(bjM.B, bjM.F,'Ts', TS, 'variable', 'z^-1');
H_bj1 = tf(bjM.C, bjM.D,'Ts', TS, 'variable', 'z^-1');
G_oe1 = tf(oeM.B, oeM.F,'Ts', TS, 'variable', 'z^-1');
H_oe1 = tf(oeM.C, oeM.D,'Ts', TS, 'variable', 'z^-1');

% SNR = 0.1
SNR = 0.1;

var_v = var(y0_det)/SNR;
e = randn(N,1);
H0z_var = lsim(H0z, e, step);
var_e = var_v / var(H0z_var);

e0 = sqrt(var_e) * randn(N,1);
y0_sim = lsim(G0z*S0z, r) + lsim(S0z*H0z, e0);
u0 = lsim(S0z, r) - lsim(R0z*S0z*H0z, e0);
data = iddata(y0_sim(500:N), u0(500:N), TS, 'Domain', 'Time');

bjM = bj(data,[2 1 1 4 3]);
oeM = oe(data, [2 4 3]);

G_bj2 = tf(bjM.B, bjM.F,'Ts', TS, 'variable', 'z^-1');
H_bj2 = tf(bjM.C, bjM.D,'Ts', TS, 'variable', 'z^-1');
G_oe2 = tf(oeM.B, oeM.F,'Ts', TS, 'variable', 'z^-1');
H_oe2 = tf(oeM.C, oeM.D,'Ts', TS, 'variable', 'z^-1');


% SNR = 0.1 and N = 100000
SNR = 0.1;
N = 100000;
step = (0:TS:(N-1)*TS);

r2 = idinput(N, 'prbs', band, range);
r = lsim(R0z, r2, step);
y0_det = lsim(Lz, r2, step);

var_v = var(y0_det)/SNR;
e = randn(N,1);
H0z_var = lsim(H0z, e, step);
var_e = var_v / var(H0z_var);

e0 = sqrt(var_e) * randn(N,1);
y0_sim = lsim(G0z*S0z, r) + lsim(S0z*H0z, e0);
u0 = lsim(S0z, r) - lsim(R0z*S0z*H0z, e0);
data = iddata(y0_sim(500:N), u0(500:N), TS, 'Domain', 'Time');

bjM = bj(data,[2 1 1 4 3])
oeM = oe(data, [2 4 3])

G_bj3 = tf(bjM.B, bjM.F,'Ts', TS, 'variable', 'z^-1');
H_bj3 = tf(bjM.C, bjM.D,'Ts', TS, 'variable', 'z^-1');
G_oe3 = tf(oeM.B, oeM.F,'Ts', TS, 'variable', 'z^-1');
H_oe3 = tf(oeM.C, oeM.D,'Ts', TS, 'variable', 'z^-1');

%% plotting
[magnitude] = bode(G0z, pulse);
mod_G0z = squeeze(magnitude);
[magnitude] = bode(H0z, pulse);
mod_H0z = squeeze(magnitude);


% First row
[magnitude] = bode(G_bj1, pulse);
mod_Gbj1 = squeeze(magnitude);
[magnitude] = bode(H_bj1, pulse);
mod_Hbj1 = squeeze(magnitude);

[magnitude] = bode(G_oe1, pulse);
mod_Goej1 = squeeze(magnitude);
[magnitude] = bode(H_oe1, pulse);
mod_Hoej1 = squeeze(magnitude);


figure
tiledlayout(3,2)

nexttile
semilogx(pulse, db(mod_G0z), 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$|G_0(z)|$');
hold on
semilogx(pulse, db(mod_Gbj1), '--', 'Color', 'g', 'LineWidth', 2, 'DisplayName', '$|G_{bj}(z)|$');
hold on
semilogx(pulse, db(mod_Goej1), ':','Color', 'b', 'LineWidth', 2, 'DisplayName', '$|G_{oe}(z)|$');
hold on
set(legend('Interpreter','Latex'))
s1 = legend
s1.Location = 'southwest';
title('Estimate of $G_0(z)$ with SNR = 100', 'Interpreter', 'latex')
ylabel("Magnitude [dB]");
xlabel("Frequency [rad/s]");
xlim([10^0, 5*10^1]);
ylim([-40, 40]);
yticks(-40:20:40);


nexttile
semilogx(pulse, db(mod_H0z), 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$|H_0(z)|$');
hold on
semilogx(pulse, db(mod_Hbj1), '--', 'Color', 'g', 'LineWidth', 2, 'DisplayName', '$|H_{bj}(z)|$');
hold on
semilogx(pulse, db(mod_Hoej1), ':','Color', 'b', 'LineWidth', 2, 'DisplayName', '$|H_{oe}(z)|$');
hold on
set(legend('Interpreter','Latex'))
s1 = legend
s1.Location = 'southwest';
title('Estimate of $H_0(z)$ with SNR = 100', 'Interpreter', 'latex')
ylabel("Magnitude [dB]");
xlabel("Frequency [rad/s]");
xlim([10^0, 5*10^1]);
ylim([-40, 40]);
yticks(-40:20:40);


% 2 row 
[magnitude] = bode(G_bj2, pulse);
mod_Gbj2 = squeeze(magnitude);
[magnitude] = bode(H_bj2, pulse);
mod_Hbj2 = squeeze(magnitude);

[magnitude] = bode(G_oe2, pulse);
mod_Goej2 = squeeze(magnitude);
[magnitude] = bode(H_oe2, pulse);
mod_Hoej2 = squeeze(magnitude);



nexttile
semilogx(pulse, db(mod_G0z), 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$|G_0(z)|$');
hold on
semilogx(pulse, db(mod_Gbj2), '--', 'Color', 'g', 'LineWidth', 2, 'DisplayName', '$|G_{bj}(z)|$');
hold on
semilogx(pulse, db(mod_Goej2), ':','Color', 'b', 'LineWidth', 2, 'DisplayName', '$|G_{oe}(z)|$');
hold on
set(legend('Interpreter','Latex'))
s1 = legend
s1.Location = 'southwest';
title('Estimate of $G_0(z)$ with SNR = 0.1', 'Interpreter', 'latex')
ylabel("Magnitude [dB]");
xlabel("Frequency [rad/s]");
xlim([10^0, 5*10^1]);
ylim([-40, 40]);
yticks(-40:20:40);


nexttile
semilogx(pulse, db(mod_H0z), 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$|H_0(z)|$');
hold on
semilogx(pulse, db(mod_Hbj2), '--', 'Color', 'g', 'LineWidth', 2, 'DisplayName', '$|H_{bj}(z)|$');
hold on
semilogx(pulse, db(mod_Hoej2), ':','Color', 'b', 'LineWidth', 2, 'DisplayName', '$|H_{oe}(z)|$');
hold on
set(legend('Interpreter','Latex'))
s1 = legend
s1.Location = 'southwest';
title('Estimate of $H_0(z)$ with SNR = 0.1', 'Interpreter', 'latex')
ylabel("Magnitude [dB]");
xlabel("Frequency [rad/s]");
xlim([10^0, 5*10^1]);
ylim([-40, 40]);
yticks(-40:20:40);

% 3 row 
[magnitude] = bode(G_bj3, pulse);
mod_Gbj3 = squeeze(magnitude);
[magnitude] = bode(H_bj3, pulse);
mod_Hbj3 = squeeze(magnitude);

[magnitude] = bode(G_oe3, pulse);
mod_Goej3 = squeeze(magnitude);
[magnitude] = bode(H_oe3, pulse);
mod_Hoej3 = squeeze(magnitude);

nexttile
semilogx(pulse, db(mod_G0z), 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$|G_0(z)|$');
hold on
semilogx(pulse, db(mod_Gbj3), '--', 'Color', 'g', 'LineWidth', 2, 'DisplayName', '$|G_{bj}(z)|$');
hold on
semilogx(pulse, db(mod_Goej3), ':','Color', 'b', 'LineWidth', 2, 'DisplayName', '$|G_{oe}(z)|$');
hold on
set(legend('Interpreter','Latex'))
s1 = legend
s1.Location = 'southwest';
title('Estimate of $G_0(z)$ with SNR = 0.1 with N = 100000', 'Interpreter', 'latex')
ylabel("Magnitude [dB]");
xlabel("Frequency [rad/s]");
xlim([10^0, 5*10^1]);
ylim([-40, 40]);
yticks(-40:20:40);


nexttile
semilogx(pulse, db(mod_H0z), 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$|H_0(z)|$');
hold on
semilogx(pulse, db(mod_Hbj3), '--', 'Color', 'g', 'LineWidth', 2, 'DisplayName', '$|H_{bj}(z)|$');
hold on
semilogx(pulse, db(mod_Hoej3), ':','Color', 'b', 'LineWidth', 2, 'DisplayName', '$|H_{oe}(z)|$');
hold on
set(legend('Interpreter','Latex'))
s1 = legend
s1.Location = 'southwest';
title('Estimate of $H_0(z)$ with SNR = 0.1 with N = 100000', 'Interpreter', 'latex')
ylabel("Magnitude [dB]");
xlabel("Frequency [rad/s]");
xlim([10^0, 5*10^1]);
ylim([-40, 40]);
yticks(-40:20:40);
