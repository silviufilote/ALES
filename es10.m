clear;
clc;
rng("default");
close all;
load("data_mimo_moesp.mat");
set(0, 'defaultAxesFontSize', 9)
set(0, 'DefaultLineLineWidth', 2);
set(0, 'defaultAxesFontSize', 9)
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultlegendInterpreter','latex')


%% Setup
N= 30;
time = [0:1:N-1];

%% Deterministic data
p=2;
mu=3;
A = [0.7 0; 0 0.3];
B = [-0.1 0.2 0.6; 0.9 -0.5 -0.4];
C = [0.5 0.2; 0.6 -0.8];
D = zeros(p, mu);

G0z = ss(A,B,C,D,1);
y_det = impulse(G0z,N-1);

y_det1 = y_det(:,:,1);
y_det2 = y_det(:,:,2);
y_det3 = y_det(:,:,3);

%% estimating MOESP with function
[ss_out,ssmat]= moesp(y',u',5);
[A,B,C,D] = ssmat(4);
G_hat = ss(A,B,C,D,1);
y_sim_det = impulse(G_hat, N-1);

y1_sim_det1= y_sim_det(:,1,1);
y1_sim_det2= y_sim_det(:,1,2);
y1_sim_det3= y_sim_det(:,1,3);

y2_sim_det1= y_sim_det(:,2,1);
y2_sim_det2= y_sim_det(:,2,2);
y2_sim_det3= y_sim_det(:,2,3);

%% Plotting
figure('Name','Noisy data')
tiledlayout(2,3)
nexttile 
plot(time, y1_sim_det1(:,1), 'b',  'DisplayName', "Estimated impulse response");
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
plot(time, y1_sim_det2(:,1), 'b');
hold on;

plot(time, y_det2(:,1), 'k--');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.07]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 2 to output 1");


nexttile
plot(time, y1_sim_det3(:,1), 'b');
hold on;
plot(time, y_det3(:,1), 'k--');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.3]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 3 to output 1");

nexttile 
plot(time, y2_sim_det1(:,1), 'b');

hold on;
plot(time, y_det1(:,2), 'k--');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([-1, 0]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 1 to output 2");

nexttile 
plot(time, y2_sim_det2(:,1), 'b');
hold on;
plot(time, y_det2(:,2), 'k--');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.6]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 1 to output 2");

nexttile 
plot(time, y2_sim_det3(:,1), 'b');
hold on;
plot(time, y_det3(:,2), 'k--');
xlabel('Time instants');
ylabel('Amplitude');
grid on;
ylim([0, 0.8]);
xlim([0,30]); xticks(0:10:time(end));
title("From input 1 to output 2");

%% Function for moesp
function [ss,ssfun]=moesp(y,u,d)

[ndat,ny]=size(y);
[mdat,nu]=size(u);
if ndat~=mdat
    error('Y and U have different length.')
end

% block Hankel matrix
N=ndat-d+1;
Y = zeros(d*ny,N);
U = zeros(d*nu,N);
% sN = sqrt(N);
sy = y';  %/sN
su = u';  %/sN
for s=1:d
    Y((s-1)*ny+1:s*ny,:)=sy(:,s:s+N-1);
    U((s-1)*nu+1:s*nu,:)=su(:,s:s+N-1);
end

% RQ decomposition
R=triu(qr([U;Y]'))';
R=R(:,1:d*(ny+nu));

% SVD
R22 = R(d*nu+1:end,d*nu+1:end);
[U1,S1]=svd(R22);

plot(svd(R22), 'b--o');  %% scala asse y? 
hold on;
xlabel('Singular value number');
ylabel('singular value');
grid on;
xlim([1,10]);
title("Singular values");

% sigular value
ss = diag(S1);

ssfun = @ssmat;

    function [A,B,C,D]=ssmat(n)
        % C and A
        Ok = U1(:,1:n)*diag(sqrt(ss(1:n)));
        C=Ok(1:ny,:);
        A=Ok(1:ny*(d-1),:)\Ok(ny+1:d*ny,:);

        % B and D
        L1 = U1(:,n+1:end)';
        R11 = R(1:d*nu,1:d*nu);
        R21 = R(d*nu+1:end,1:d*nu);
        M1 = L1*R21/R11;
        m = ny*d-n;
        M = zeros(m*d,nu);
        L = zeros(m*d,ny+n);
        for k=1:d
            M((k-1)*m+1:k*m,:)=M1(:,(k-1)*nu+1:k*nu);
            L((k-1)*m+1:k*m,:)=[L1(:,(k-1)*ny+1:k*ny) L1(:,k*ny+1:end)*Ok(1:end-k*ny,:)];
        end
        DB=L\M;
        D=DB(1:ny,:);
        B=DB(ny+1:end,:);
    end
end


