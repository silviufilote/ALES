clear;
clc;
rng("default");
close all;

%% setup
N = 10000;                          % number of data collected
TS = 0.05;                          % sample time [s]
fs = 1/TS;                          % sampling frequency [Hz]
pulse = 1e0:(1e-2):(fs*pi);         % 1:0.01:
SNR = 1e8;                          % Signal noise ratio

% creates special variable z that you can use in a rational expression 
% to create a discrete-time transfer function model. TS = Sample time.
G0z = tf([0 0 0 0.28261 0.50666],[1 -1.41833 1.58939 -1.31608 0.88642],'Ts', TS, 'variable', 'z^-1');
R0z = tf([0.32905 -0.59771 0.70728 -0.64010 0.46499 -0.11769], [1 -1], 'Ts', TS, 'variable', 'z^-1');
S0z = 1/(1 + G0z*R0z);

%% RPBS
range = [-2, 2];
band = [0 0.1];                    

% T (periodo) = N che sono i dati che voglio collezionare
% 1 = Number of input channels in generated signal, specified as a real positive integer.
% 1 = Number of periods in generated signal.
r2 = idinput([N 1 1], 'prbs', band, range);
r = lsim(R0z, r2);

%% Non parametric estimation 
step = (0:TS:(N-1)*TS);
Fz = feedback(R0z*G0z, 1);

% Plots the simulated time response of the dynamic system model Gz
% Gz = Dynamic system, specified as a SISO or MIMO 
% r2 = Input signal for simulation, specified as a vector for single-input systems, and an array for multi-input systems.
y0 = lsim(Fz, r2, step);

% lambda noise
var_v = var(y0)/SNR;
e1 = sqrt(var_v) * randn(N,1);
y1 = lsim(G0z*S0z, r) + lsim(S0z, e1);
u1 = lsim(S0z, r) - lsim(R0z*S0z, e1);

% signals
t1 = 6/TS;
t2 = 20/TS;

% Encapsulate input and output measurement data for the system you want to
% identify Domain = "Time" => Data is in the time domain
data = iddata(y1(500:N), u1(500:N), TS, 'Domain', 'Time');

% estimates the input-to-output frequency response G(ω) and noise spectrum Φυ of the general linear model
G_hat = spafdr(data, [], pulse);


%% Plotting Estimates G0(z), 1/R(z), G_hat(z)

% The plot displays the magnitude (in dB) and phase (in degrees) of the system response as a function of frequency
% [mag,phase] = bode(sys,w) returns the response data at the frequencies specified by w.
% If w is a cell array of the form {wmin,wmax}, then wout contains frequencies ranging between wmin and wmax.
% If w is a vector of frequencies, then wout = w.

[magnitude] = bode(G0z, pulse);
mod_G0z = squeeze(magnitude);

[magnitude] = bode(R0z, pulse);
mod_R0z = squeeze(magnitude);

[magnitude] = bode(G_hat, pulse);
mod_G_hat_1 = squeeze(magnitude);

SNR = 1e-8;                          %%%%  Signal noise ratio
var_v = var(y0)/SNR;
e2 = sqrt(var_v) * randn(N,1);
y2 = lsim(G0z*S0z, r) + lsim(S0z, e2);
u2 = lsim(S0z, r) - lsim(R0z*S0z, e2);
data = iddata(y2(500:N), u2(500:N), TS, 'Domain', 'Time');
G_hat = spafdr(data, [], pulse);

[magnitude] = bode(G_hat, pulse);
mod_G_hat_2 = squeeze(magnitude);

SNR = 5;                            %%%%  Signal noise ratio
var_v = var(y0)/SNR;
e3 = sqrt(var_v) * randn(N,1);
y3 = lsim(G0z*S0z, r) + lsim(S0z, e3);
u3 = lsim(S0z, r) - lsim(R0z*S0z, e3);
data = iddata(y3(500:N), u3(500:N), TS, 'Domain', 'Time');
G_hat = spafdr(data, []);

[magnitude, phase] = bode(G_hat, pulse);
mod_G_hat_3 = squeeze(magnitude);

%%  first row

figure
tiledlayout(3,2)

nexttile
semilogx(pulse, db(-1./mod_R0z), 'b', 'LineWidth', 2, 'DisplayName', '$|-1/R(z)|$');
hold on 
semilogx(pulse, db(mod_G0z), 'k--', 'LineWidth', 2, 'DisplayName', '$|G_0(z)|$');
hold on 
semilogx(pulse, db(mod_G_hat_1), ':', 'Color', 'g', 'LineWidth', 2, 'DisplayName', '$|\hat{G}(e^{j\omega})|$');
set(legend('Interpreter','Latex'))
s1 = legend
s1.Location = 'southwest';

title('Estimate of $G_0(z)$ with SNR = 10^8', 'Interpreter', 'latex')
ylabel("Magnitude [dB]");
xlabel("Frequency [rad/s]");
ylim([-40, 40]);
yticks(-40:20:40);
xlim([10^0, 5*10^1]);


nexttile
hold on
plot([6:TS:20], r2(t1:t2), "color", "black", "linewidth", 2, "linestyle", "--")
plot([6:TS:20], y1(t1:t2), "color", "blue", "linewidth", 2)
plot([6:TS:20], u1(t1:t2), "color", "red", "linewidth", 2)


%% second row

nexttile
semilogx(pulse, db(-1./mod_R0z), 'b', 'LineWidth', 2, 'DisplayName', '$|-1/R(z)|$');
hold on 
semilogx(pulse, db(mod_G0z), 'k--', 'LineWidth', 2, 'DisplayName', '$|G_0(z)|$');
hold on 
semilogx(pulse, db(mod_G_hat_2),':', 'Color', 'g', 'LineWidth', 2, 'DisplayName', '$|\hat{G}(e^{j\omega})|$');
set(legend('Interpreter','Latex'))
s1 = legend
s1.Location = 'southwest';

title('Estimate of $G_0(z)$ with SNR = 10^-8', 'Interpreter', 'latex')
ylabel("Magnitude [dB]");
xlabel("Frequency [rad/s]");
ylim([-40, 40]);
yticks(-40:20:40);
xlim([10^0, 5*10^1]);

nexttile
hold on
plot([6:TS:20], r2(t1:t2), "color", "black", "linewidth", 2, "linestyle", "--")
plot([6:TS:20], y2(t1:t2), "color", "blue", "linewidth", 2)
plot([6:TS:20], u2(t1:t2), "color", "red", "linewidth", 2)



%% third row
nexttile
semilogx(pulse, db(-1./mod_R0z), 'b', 'LineWidth', 2, 'DisplayName', '$|-1/R(z)|$');
hold on 
semilogx(pulse, db(mod_G0z), 'k--', 'LineWidth', 2, 'DisplayName', '$|G_0(z)|$');
hold on 
semilogx(pulse, db(mod_G_hat_3), ':', 'Color', 'g', 'LineWidth', 2, 'DisplayName', '$|\hat{G}(e^{j\omega})|$');
set(legend('Interpreter','Latex'))
s1 = legend
s1.Location = 'southwest';

title('Estimate of $G_0(z)$ with SNR = 5', 'Interpreter', 'latex')
ylabel("Magnitude [dB]");
xlabel("Frequency [rad/s]");
ylim([-40, 40]);
yticks(-40:20:40);
xlim([10^0, 5*10^1]);


nexttile
hold on
plot([6:TS:20], r2(t1:t2), "color", "black", "linewidth", 2, "linestyle", "--")
plot([6:TS:20], y3(t1:t2), "color", "blue", "linewidth", 2)
plot([6:TS:20], u3(t1:t2), "color", "red", "linewidth", 2)
