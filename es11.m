clearvars -except u_arrow 
close all
rng('default');
%% data generation

% sampling time and number of samples
Ts = 1;
N = 200;

% input, disturbances and fault orders
mu = 2;
md = 3;
mf = 1;

% output and state order
p = 3;
n = 4;

% input, disturbances and fault generation
u = mvnrnd(zeros(mu,1), eye(mu),N);
d = mvnrnd(zeros(md,1), eye(md),N);
f = [zeros(N/2,mf); ones(N/2,mf)];

%% system declaration

% system matrixes
A = [0.5 -0.7 0.7 0; 0 0.8 0.06 0; -1 0 0 0.1; 0 0 -0.1 0.4];
B = [0 0; 1 0; 0 1; 0 0];
C = [0 0 1 0; 0 1 0 0; 0 0 0 1]; 
D = zeros(p, mu);

% disturbances matrixes
Bd = [0 1 0; 2 0 0; 0 0 0; 0 0 0];
Dd = zeros(p, md);

% faults matrixes
Bf = zeros(n, mf);
Df = [1; 0; 0];

% real system
sys = ss(A, B, C, D, Ts);
[y, t, x] = lsim(sys, u);

% system with distrubances 
B_d = [B Bd];
D_d = [D Dd];
u_d = [u d];

sys_d = ss(A, B_d, C, D_d, Ts);
[y_d, t, x_d] = lsim(sys_d, u_d);

% system with disurbances and faults
B_df = [B_d Bf];
D_df = [D_d Df];
u_df = [u_d f];

sys_df = ss(A, B_df, C, D_df, Ts);
y_df = lsim(sys_df, u_df);

%% plot of state and output (real, disturbed and faulty)

% state plot
figure;
for i = 1:n
    subplot(n/2,n/2,i);
    hold on;
    plot(t, x_d(:,i), 'linewidth', 2, 'color', 'blue');
    plot(t, x(:,i), 'linewidth', 2, 'linestyle', ':', 'color', 'black');
    text = ['State x_{', num2str(i), '}(t)'];
    title(text, 'Interpreter', 'Tex');
    legend('Disturbed', 'True');
    xlabel('Time [s]');
    ylim([-10,10])
end

% output plot
figure;
for i = 1:p
    subplot(2,p,i);
    hold on;
    plot(t, y_d(:,i), 'linewidth', 2, 'color', 'blue');
    plot(t, y(:,i), 'linewidth', 2, 'linestyle', ':', 'color', 'black');
    text = ['Output y_{', num2str(i), '}(t)'];
    title(text, 'Interpreter', 'Tex');
    legend('Disturbed', 'True');
    xlabel('Time [s]');
    ylim([-10,10])
end
for i = 1:p
    subplot(2,p,i+p);
    hold on;
    plot(t, y_df(:,i), 'linewidth', 2, 'color', 'red');
    text = ['Output y_{', num2str(i), '}(t) faulty'];
    title(text, 'Interpreter', 'Tex');
    xlabel('Time [s]');
    ylim([-10,10])
end

%% residual generator design
s = 1;

% extended observability matrix
O = zeros(p*(s+1), n);
O(1:p,:) = C;
for i = 1:s
    O(i*p+1:(i+1)*p,:) = O((i-1)*p+1:i*p,:)*A;
end

% extended disturbances matrix Ds+1
Ds = zeros(p*(s+1), md*(s+1));
Ds(1:p, 1:md) = Dd;
for i = 1:s
    Ds(i*p+1:(i+1)*p,1:md) = C*(A^(i-1))*Bd;
end
for i = 1:s
    Ds(i*p+1:end,i*md+1:(i+1)*md) = Ds(1:size(Ds,1)-i*p,1:md);
end
% vector from the parity space
V = null([O Ds]');
v = V(:,1)';

% vector from the parity space without decoupling
V_nd = null(O');
v_nd = V_nd(:,1)';

% T matrix
T = zeros(p*(s+1), mu*(s+1));
T(1:p, 1:mu) = D;
for i = 1:s
    T(i*p+1:(i+1)*p,1:mu) = C*(A^(i-1))*B;
end
for i = 1:s
    T(i*p+1:end,i*mu+1:(i+1)*mu) = T(1:size(T,1)-i*p,1:mu);
end

%input and output reshaping
Yr = reshape(y_df',[],1);
Ur = reshape(u',[],1);

%computing 'arrow' matrixes
Ys = zeros(p*(s+1), N-s);
Us = zeros(mu*(s+1), N-s);

for i = 1:(N - s)
    Ys(:,i) = Yr(p*(i-1)+1:p*(i-1+s+1));
    Us(:,i) = Ur(mu*(i-1)+1:mu*(i-1+s+1));
end

%residual signal
r = abs(v*(Ys - T*Us));
r_nd = abs(v_nd*(Ys - T*Us));

% plot of residuals signals
figure;
subplot(2,1,1);
plot(t(s+1:end), r_nd, 'linewidth', 2, 'color', 'black');
ylabel('|r(t)|', 'Interpreter', 'Tex')
xlabel('Time [s]')
title('Residual not decoupled')

subplot(2,1,2);
plot(t(s+1:end), r, 'linewidth', 2, 'color', 'black');
ylabel('|r(t)|', 'Interpreter', 'Tex')
xlabel('Time [s]')
title('Residual decoupled')