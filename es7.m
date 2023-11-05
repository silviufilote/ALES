clear;
clc;
rng("default");
close all;

set(0, 'defaultAxesFontSize', 9)
set(0, 'DefaultLineLineWidth', 2);
set(0, 'defaultAxesFontSize', 9)
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultlegendInterpreter','latex')


%% Function run
projectionFunc(10000, 1)
projectionFunc(10000, 2)

function projectionFunc(N, caseFigure)
    %% Setup
    SNR = 0.1;
    TS = 0.05;
    fs = 1/TS; 
    pulse = 1e0:(1e-2):(fs*pi); 
    
    % where z is the variable in the Laplace domain
    G0z = tf([0 0 0 0.28261 0.50666],[1 -1.41833 1.58939 -1.31608 0.88642],'Ts', TS, 'variable', 'z^-1');
    R0z = tf([0.32905 -0.59771 0.70728 -0.64010 0.46499 -0.11769], [1 -1], 'Ts', TS, 'variable', 'z^-1');
    H0z = tf([1 0.134], [1 -0.789], 'Ts', TS, 'variable', 'z^-1');
    S0z = 1/(1 + G0z*R0z);

    %% RPBS
    Range = [-2, 2];
    Band = [0 1]; 
    step = (0:TS:(N-1)*TS);
    input_channels = 1;
    nperiods = 1;
    
    r2 = idinput([N input_channels nperiods], 'prbs', Band, Range);
    r = lsim(R0z, r2, step);

    %% Simulate true signal 
    Fz = feedback(R0z*G0z, 1);
    y0_det = lsim(Fz, r2, step);
    
    e = randn(N,1);
    v = lsim(H0z, e, step);
    var_v = var(v);
    varnoise = var(y0_det)/(SNR * var_v);
    e0 = sqrt(varnoise) * randn(N, 1);
    
    y0_sim = lsim(G0z*S0z, r) + lsim(S0z*H0z, e0);
    u0 = lsim(S0z, r) - lsim(R0z*S0z*H0z, e0);
    data = iddata(u0(500:end), r(500:end), TS);

    %% direct identification method -> PEM estimation (open loop)
    
    nb = 5;
    nc = 10;
    nd = 10;
    nf = 9;
    nk = 0;
    
    bjM = bj(data, [nb+1 nc nd nf nk]);
    S0_hat_bj = tf(bjM.B, bjM.F,'Ts', TS, 'variable', 'z^-1');
    
    u_r = lsim(S0_hat_bj, r, step);
    data = iddata(y0_sim(500:end), u_r(500:end), TS);
    
    nb = 1;
    nc = 6;
    nd = 10;
    nf = 4;
    nk = 3;
    
    bjM = bj(data, [nb+1 nc nd nf nk]);
    G_bj1 = tf(bjM.B, bjM.F,'Ts', TS, 'variable', 'z^-1');
    W_bj1 = tf(bjM.C, bjM.D,'Ts', TS, 'variable', 'z^-1');
    H_bj1 = minreal(W_bj1 / S0_hat_bj);
    
    % oe
    
    oeM = oe(data, [nb+1 nf nk]);
    G_oe1 = tf(oeM.B, oeM.F,'Ts', TS, 'variable', 'z^-1');
    W_oe1 = tf(1, TS, 'variable', 'z^-1');
    H_oe1 = minreal(W_oe1 / S0_hat_bj);

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
    
    switch caseFigure
        case 1
            figure
            tiledlayout(2,2)
            
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
            title('Estimate of $G_0(z)$ with SNR = 0.1 with N = 10000', 'Interpreter', 'latex')
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
            title('Estimate of $H_0(z)$ with SNR = 0.1 with N = 10000', 'Interpreter', 'latex')
            ylabel("Magnitude [dB]");
            xlabel("Frequency [rad/s]");
            xlim([10^0, 5*10^1]);
            ylim([-40, 40]);
            yticks(-40:20:40);

       case 2
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
            title('Estimate of $G_0(z)$ with SNR = 0.1 with N = 100000', 'Interpreter', 'latex')
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
            title('Estimate of $H_0(z)$ with SNR = 0.1 with N = 100000', 'Interpreter', 'latex')
            ylabel("Magnitude [dB]");
            xlabel("Frequency [rad/s]");
            xlim([10^0, 5*10^1]);
            ylim([-40, 40]);
            yticks(-40:20:40);

       otherwise
            disp('error')
    end
end
