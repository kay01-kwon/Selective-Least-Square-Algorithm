close all
clear all
load simulation2.mat

params = gyroparams;

dt = 0.01;

dzdt = s_RK4(3,:);
z = s_RK4(6,:);

acc = zeros(length(s_RK4(7,:)),3);

imu = imuSensor('SampleRate', 100, 'Gyroscope', params);
imu.Gyroscope.NoiseDensity = deg2rad(0.0035);
imu.Gyroscope.Resolution = deg2rad(0.00625);
[~, angular_Vel] = imu(acc, s_RK4(7:9,:)');

p = angular_Vel(:,1)';
q = angular_Vel(:,2)';
r = angular_Vel(:,3)';

% p = s_RK4(7,:);
% q = s_RK4(8,:);
% r = s_RK4(9,:);

T = u_vec(1,:)';
M_x = u_vec(2,:)';
M_y = u_vec(3,:)';

N = length(z);

% Estimation Parameter Setup
J_xx_hat = [];
J_xx_hat_vec = [];
J_xx_var = 0.005^2;
sigma1 = 0.001;
sigma2 = 0.0001;
J_xx_var_vec = [J_xx_var];

y_off_hat = [];
y_off_hat_vec = [];
y_off_var = 0.010^2;
y_off_var_vec = [y_off_var];

J_yy_hat = [];
J_yy_var = 0.005^2;
J_yy_var_vec = [J_yy_var];

x_off_hat = [];
x_off_hat_vec = [];
x_off_var = 0.010^2;
x_off_var_vec = [x_off_var];

x = [p(2);-0.0981];
dpdt_ = 0;
dqdt_ = 0;

b_xx_vec = [];
b_yy_vec = [];

% Initial Data sampling
for i = 2:4
    dpdt(i) = (p(i+1) - p(i))/dt;
    dqdt(i) = (q(i+1) - q(i))/dt;

    A_xx(i-1,:) = [dpdt(i) p(i)*r(i) T(i)];
    A_yy(i-1,:) = [dqdt(i) q(i)*r(i) -T(i)];    
end

W_xx = A_xx'*A_xx;
W_yy = A_yy'*A_yy;

b_xx = A_xx'*M_x(2:4);
b_yy = A_yy'*M_y(2:4);

b_xx_vec = [b_xx_vec;b_xx];
b_yy_vec = [b_yy_vec;b_yy];

L_xx_vec = [];
L_yy_vec = [];

time_stamp_J_xx = [];
time_stamp_J_yy = [];


for i = 4:N-1
    dpdt(i) = (p(i+1) - p(i))/dt;
    dqdt(i) = (q(i+1) - q(i))/dt;

    a_xx = [dpdt(i) p(i)*r(i) T(i)];
    a_yy = [dqdt(i) q(i)*r(i) -T(i)];
    
    [W_xx b_xx] = FimUpdate(W_xx,b_xx,a_xx,M_x(i-2));
    [W_yy b_yy] = FimUpdate(W_yy,b_yy,a_yy,M_y(i-2));
    
    L_xx = inv(W_xx)*b_xx;
    L_yy = inv(W_yy)*b_yy;
    
    L_xx_vec = [L_xx_vec L_xx];
    L_yy_vec = [L_yy_vec L_yy];
    
    [V_xx,D_xx] = eig(W_xx);
    [V_xx_obs,D_xx_obs] = CheckSVD(V_xx,D_xx);
    obs_xx_index = CheckObsParam(V_xx_obs);
    solution_xx(i-3,:) = GetRecursiveSolution(V_xx_obs,D_xx_obs,b_xx);
    
    [V_yy,D_yy] = eig(W_yy);
    [V_yy_obs,D_yy_obs] = CheckSVD(V_yy,D_yy);
    obs_yy_index = CheckObsParam(V_yy_obs);
    solution_yy(i-3,:) = GetRecursiveSolution(V_yy_obs,D_yy_obs,b_yy);

    n_xx = length(obs_xx_index);
    n_yy = length(obs_yy_index);    

    % Estimation along x axis and update
    for k = 1:n_xx
        if obs_xx_index(k) == 1
            if isempty(J_xx_hat) == 1
                J_xx_hat = solution_xx(i-3,1);
                J_xx_hat_vec = [J_xx_hat];
            else
                K_J_xx = J_xx_var/(J_xx_var + sigma1);
                J_xx_hat = J_xx_hat + K_J_xx*(solution_xx(i-3,1) - J_xx_hat);
                J_xx_var = (1-K_J_xx)*J_xx_var;
                J_xx_hat_vec = [J_xx_hat_vec J_xx_hat];
                J_xx_var_vec = [J_xx_var_vec J_xx_var];
                time_stamp_J_xx = [time_stamp_J_xx t(i)];
            end
        elseif obs_xx_index(k) == 3
            if isempty(y_off_hat) == 1
                y_off_hat = solution_xx(i-3,3);
            else
                K_y_off = y_off_var/(y_off_var + sigma2);
                y_off_hat = y_off_hat + K_y_off*(solution_xx(i-3,3) - y_off_hat);
                y_off_var = (1-K_y_off)*y_off_var;
                y_off_hat_vec = [y_off_hat_vec y_off_hat];
                y_off_var_vec = [y_off_var_vec y_off_var];
            end
        end
    end
    % End of For statement

    % Estimation along y axis and update
    for k = 1:n_yy
        if obs_yy_index(k) == 1
            if isempty(J_yy_hat) == 1
                J_yy_hat = solution_yy(i-3,1);
                J_yy_hat_vec = [J_yy_hat];
            else
                K_J_yy = J_yy_var/(J_yy_var + sigma1);
                J_yy_hat = J_yy_hat + K_J_yy*(solution_yy(i-3,1) - J_yy_hat);
                J_yy_var = (1-K_J_yy)*J_yy_var;
                J_yy_hat_vec = [J_yy_hat_vec J_yy_hat];
                J_yy_var_vec = [J_yy_var_vec J_yy_var];
                time_stamp_J_yy = [time_stamp_J_yy t(i)];
            end
        elseif obs_yy_index(k) == 3
            if isempty(x_off_hat) == 1
                x_off_hat = solution_yy(i-3,3);
            else
                K_x_off = x_off_var/(x_off_var + sigma2);
                x_off_hat = x_off_hat + K_x_off*(solution_yy(i-3,3) - x_off_hat);
                x_off_var = (1-K_x_off)*x_off_var;
                x_off_hat_vec = [x_off_hat_vec x_off_hat];
                x_off_var_vec = [x_off_var_vec x_off_var];
            end
        end
    end
    
    b_xx_vec = [b_xx_vec;b_xx];
    b_yy_vec = [b_yy_vec;b_yy];
end

figure(1)
subplot(2,2,1)
plot(t(2:end-1),dpdt(2:end))
hold on
plot(t(3:end),dsdt(7,3:end),'--','LineWidth',2)
title('dpdt - t ')
legend('dpdt(diff)','dpdt(ODE)')
xlabel('time (s)')
ylabel('dpdt (rad/s^2)')
grid on

subplot(2,2,2)
plot(dqdt(2:end))
hold on
plot(dsdt(8,3:end),'--','LineWidth',2)
title('dpdt - t ')
legend('dqdt(diff)','dqdt(ODE)')
xlabel('time (s)')
ylabel('dqdt (rad/s^2)')
grid on

subplot(2,2,3)
plot(t(2:end-1),p(3:end))
hold on;
plot(t(3:end),s_RK4(7,3:end),'--')
title('p - t')
legend('p(IMU)','p(ODE)')
xlabel('time (s)')
ylabel('p (rad/s)')
grid on

subplot(2,2,4)
plot(t(3:end),q(3:end))
hold on;
plot(t(3:end),s_RK4(8,3:end),'--')
title('q - t ')
legend('q(IMU)','q(ODE)')
xlabel('time (s)')
ylabel('q (rad/s)')
grid on

figure()
if isempty(J_xx_hat_vec) == 0
subplot(2,2,1)
t_J_xx = 1:length(J_xx_hat_vec);
J_xx_hat_vec_plus = [t_J_xx' (J_xx_hat_vec + sqrt(J_xx_var_vec))'];
J_xx_hat_vec_minus = [t_J_xx' (J_xx_hat_vec - sqrt(J_xx_var_vec))'];
J_xx_between = [J_xx_hat_vec_minus;flipud(J_xx_hat_vec_plus)];
fill(J_xx_between(:,1),J_xx_between(:,2),[0 1 0.9])
hold on;
plot(J_xx_hat_vec,'r','LineWidth',2);
grid on;
title('$\hat{J}_{xx}-Estimation\,step$','Interpreter','latex')
xlabel('Estimation step(n)')
ylabel('$\hat{J}_{xx} (kg \cdot m^2) $','Interpreter','latex')
% str = '1.5e-3 kg m^2';
% text(t_J_xx(end),J_xx_hat_vec(end),str)
% axis([0 t_J_xx(end) 0 0.0015*1.2])

% subplot(2,2,2)
% plot(sqrt(J_xx_var_vec));
% axis([0 length(J_xx_hat_vec) 0 0.02])
% grid on;
% title('$\sigma_{\hat{J}_{xx}} (kg \cdot m^2)-Time$','Interpreter','latex')
% xlabel('Time step(n)')
% ylabel('$\sigma_{\hat{J}_{xx}} (kg \cdot m^2)$','Interpreter','latex')
end
subplot(2,2,2)
t_J_yy = 1:length(J_yy_hat_vec);
J_yy_hat_vec_plus = [t_J_yy' (J_yy_hat_vec + sqrt(J_yy_var_vec))'];
J_yy_hat_vec_minus = [t_J_yy' (J_yy_hat_vec - sqrt(J_yy_var_vec))'];
J_yy_between = [J_yy_hat_vec_minus;flipud(J_yy_hat_vec_plus)];
fill(J_yy_between(:,1),J_yy_between(:,2),[0 1 0.9])
hold on;
plot(J_yy_hat_vec,'r','LineWidth',2);
grid on;
title('$\hat{J}_{yy}-Estimation\,step$','Interpreter','latex')
xlabel('Estimation step(n)')
ylabel('$\hat{J}_{yy} (kg \cdot m^2) $','Interpreter','latex')
% str = '1.5e-3 kg m^2';
% text(t_J_yy(end),J_yy_hat_vec(end),str)


subplot(2,2,3)
t_y_off = 1:length(y_off_hat_vec);
y_off_hat_vec_plus = [t_y_off' (y_off_hat_vec + sqrt(y_off_var_vec(1:end-1)))'];
y_off_hat_vec_minus = [t_y_off' (y_off_hat_vec - sqrt(y_off_var_vec(1:end-1)))'];
y_off_between = [y_off_hat_vec_minus;flipud(y_off_hat_vec_plus)];
fill(y_off_between(:,1),y_off_between(:,2),[0 1 1]);

hold on;
plot(y_off_hat_vec,'b','LineWidth',2);
% axis([0 length(y_off_hat_vec) 0 0.02])
grid on;
title('$\hat{y}_{CG \rightarrow COM}-Estimation\,step$','Interpreter','latex')
xlabel('Estimation step(n)')
ylabel('$\hat{y}_{CG \rightarrow COM} (m)$','Interpreter','latex')

% subplot(2,2,4)
% plot(sqrt(y_off_var_vec));
% axis([0 length(y_off_var_vec) 0 0.01])
% grid on;
% title('$\sigma_{\hat{y}_{CG \rightarrow COM}(m)}-Time$','Interpreter','latex')
% xlabel('Time step(n)')
% ylabel('$\sigma_{\hat{y}_{CG \rightarrow COM}(m)}$','Interpreter','latex')


subplot(2,2,4)
t_x_off = 1:length(x_off_hat_vec);
x_off_hat_vec_plus = [t_x_off' (x_off_hat_vec + sqrt(x_off_var_vec(1:end-1)))'];
x_off_hat_vec_minus = [t_x_off' (x_off_hat_vec - sqrt(x_off_var_vec(1:end-1)))'];
x_off_between = [x_off_hat_vec_minus;flipud(x_off_hat_vec_plus)];
fill(x_off_between(:,1),x_off_between(:,2),[0 1 1]);

hold on;
plot(x_off_hat_vec,'b','LineWidth',2);
% axis([0 length(x_off_hat_vec) min(x_off_hat_vec_minus)*0.8 max(x_off_hat_vec_minus)*1.2])
grid on;
title('$\hat{x}_{CG \rightarrow COM}-Estimation\,step$','Interpreter','latex')
xlabel('Estimation step(n)')
ylabel('$\hat{x}_{CG \rightarrow COM} (m)$','Interpreter','latex')

if isempty(J_xx_hat_vec) == 0
    J_xx_hat_vec(end)
end

J_yy_hat_vec(end)
y_off_hat_vec(end)
x_off_hat_vec(end)

LSxx = inv(W_xx)*b_xx
LSyy = inv(W_yy)*b_yy


% Simple Least Square Method
figure()
subplot(2,2,1)
plot(t(1:end-4),L_xx_vec(1,:))
title('$\hat{J}_{xx}-t$','Interpreter','latex')
xlabel('t (s)')
ylabel('$\hat{J}_{xx} (kg \cdot m^2)$','Interpreter','latex')
grid on;

subplot(2,2,2)
plot(t(1:end-4),L_yy_vec(1,:))
title('$\hat{J}_{yy}-t$','Interpreter','latex')
xlabel('t (s)')
ylabel('$\hat{J}_{yy} (kg \cdot m^2)$','Interpreter','latex')
grid on

subplot(2,2,3)
plot(t(1:end-4),L_xx_vec(3,:))
title('$\hat{y}_{CG \rightarrow COM}-t$','Interpreter','latex')
xlabel('t (s)')
ylabel('$\hat{y}_{CG \rightarrow COM}(m)$','Interpreter','latex')
grid on;


subplot(2,2,4)
plot(t(1:end-4),L_yy_vec(3,:))
title('$\hat{x}_{CG \rightarrow COM} (m)-t$','Interpreter','latex')
xlabel('t (s)')
ylabel('$\hat{x}_{CG \rightarrow COM}(m)$','Interpreter','latex')

grid on;



% A = W(49:51,:);
% b = M_x(49:51);
% 
% [V D] = eig(A'*A);
% 
% [V_obs D_obs] = CheckSVD(V,D);
% obs_index = CheckObsParam(V_obs);
% solution2 = GetSolution(V_obs,D_obs,A,b)
