close all
clear all
load simulation.mat

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


A = [1 dt;
    0 1];
Q = [0.01 0;
    0 10];
H = [1 0];
R = 10;
P = 5*eye(2);
b_xx_vec = [];
b_yy_vec = [];


for i = 2:N-1
    dpdt(i) = (p(i+1) - p(i))/dt;
    dqdt(i) = (q(i+1) - q(i))/dt;
    
    
%     if i >= 2
%         p(i) = 0.3*p(i) + 0.7*p(i-1);
%         q(i) = 0.3*q(i) + 0.7*q(i-1);
%         dpdt(i) = 0.3*dpdt(i) + 0.7*dpdt(i-1);
%         dqdt(i) = 0.3*dqdt(i) + 0.7*dqdt(i-1);
%     end
    
    dpdt_ = dpdt(i);
    dqdt_ = dqdt(i);

    W_xx(i,:) = [dpdt_ p(i)*r(i) T(i)];
    W_yy(i,:) = [dqdt_ q(i)*r(i) -T(i)];
%     W_(i,:) = [dsdt(7,i+1) p(i)*r(i) T(i)];
    b_xx_vec = [b_xx_vec;M_x(i)];
    b_yy_vec = [b_yy_vec;M_y(i)];
    
    % Data sampling
    if i <= N-1 && i >= 4
        A_xx = W_xx(i-2:i,:);
        b_xx = M_x(i-2:i);
        [V_xx,D_xx] = eig(A_xx'*A_xx);
        [V_xx_obs,D_xx_obs] = CheckSVD(V_xx,D_xx);
        obs_xx_index = CheckObsParam(V_xx_obs);
        solution_xx(i-2,:) = GetSolution(V_xx_obs,D_xx_obs,A_xx,b_xx);
        
        n_xx = length(obs_xx_index);

        A_yy = W_yy(i-2:i,:);
        b_yy = M_y(i-2:i);
        [V_yy,D_yy] = eig(A_yy'*A_yy);
        [V_yy_obs,D_yy_obs] = CheckSVD(V_yy,D_yy);
        obs_yy_index = CheckObsParam(V_yy_obs);
        solution_yy(i-2,:) = GetSolution(V_yy_obs,D_yy_obs,A_yy,b_yy);
        
        n_xx = length(obs_xx_index);
        n_yy = length(obs_yy_index);
        
        % Estimation along x axis and update
        for k = 1:n_xx
            if obs_xx_index(k) == 1
                if isempty(J_xx_hat) == 1
                    J_xx_hat = abs(solution_xx(i-2,1));
                    J_xx_hat_vec = [abs(solution_xx(i-2,1))];
                else
                    K_J_xx = J_xx_var/(J_xx_var + sigma1);
                    J_xx_hat = J_xx_hat + K_J_xx*(abs(solution_xx(i-2,1)) - J_xx_hat);
                    J_xx_var = (1-K_J_xx)*J_xx_var;
                    J_xx_hat_vec = [J_xx_hat_vec J_xx_hat];
                    J_xx_var_vec = [J_xx_var_vec J_xx_var];
                end
            elseif obs_xx_index(k) == 3
                if isempty(y_off_hat) == 1
                    y_off_hat = solution_xx(i-2,3);
                else
                    K_y_off = y_off_var/(y_off_var + sigma2);
                    y_off_hat = y_off_hat + K_y_off*(solution_xx(i-2,3) - y_off_hat);
                    y_off_var = (1-K_y_off)*y_off_var;
                    y_off_hat_vec = [y_off_hat_vec y_off_hat];
                    y_off_var_vec = [y_off_var_vec y_off_var];
                end
            end
        end
        % End of For statement
        
                % Estimation along x axis and update
        for k = 1:n_yy
            if obs_yy_index(k) == 1
                if isempty(J_yy_hat) == 1
                    J_yy_hat = abs(solution_yy(i-2,1));
                    J_yy_hat_vec = [abs(solution_yy(i-2,1))];
                else
                    K_J_yy = J_yy_var/(J_yy_var + sigma1);
                    J_yy_hat = J_yy_hat + K_J_yy*(abs(solution_yy(i-2,1)) - J_yy_hat);
                    J_yy_var = (1-K_J_yy)*J_yy_var;
                    J_yy_hat_vec = [J_yy_hat_vec J_yy_hat];
                    J_yy_var_vec = [J_yy_var_vec J_yy_var];
                end
            elseif obs_yy_index(k) == 3
                if isempty(x_off_hat) == 1
                    x_off_hat = solution_yy(i-2,3);
                else
                    K_x_off = x_off_var/(x_off_var + sigma2);
                    x_off_hat = x_off_hat + K_x_off*(solution_yy(i-2,3) - x_off_hat);
                    x_off_var = (1-K_x_off)*x_off_var;
                    x_off_hat_vec = [x_off_hat_vec x_off_hat];
                    x_off_var_vec = [x_off_var_vec x_off_var];
                end
            end
        end
        
    end
    % End of Data sampling
end

figure(1)
subplot(2,3,1)
plot(dpdt(2:end))
hold on
plot(dsdt(7,3:end),'--','LineWidth',2)
legend('dpdt(diff)','dpdt(ODE)')

subplot(2,3,2)
plot(dqdt(2:end))
hold on
plot(dsdt(8,3:end),'--','LineWidth',2)
legend('dqdt(diff)','dqdt(ODE)')


subplot(2,3,4)
plot(p(3:end))
hold on;
plot(s_RK4(7,3:end),'--')

subplot(2,3,5)
plot(q(3:200))
hold on;
plot(s_RK4(8,3:200),'--')


figure()
if isempty(J_xx_hat_vec) == 0
subplot(2,2,1)
t_J_xx = 1:length(J_xx_hat_vec);
J_xx_hat_vec_plus = [t_J_xx' (J_xx_hat_vec + sqrt(J_xx_var_vec))'];
J_xx_hat_vec_minus = [t_J_xx' (J_xx_hat_vec - sqrt(J_xx_var_vec))'];
J_xx_between = [J_xx_hat_vec_minus;flipud(J_xx_hat_vec_plus)];
fill(J_xx_between(:,1),J_xx_between(:,2),[0 1 0.5])
hold on;
plot(J_xx_hat_vec,'r','LineWidth',2);
grid on;
title('$\hat{J}_{xx} (kg \cdot m^2)-Estimation\,step$','Interpreter','latex')
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
fill(J_yy_between(:,1),J_yy_between(:,2),[0 1 0.5])
hold on;
plot(J_yy_hat_vec,'r','LineWidth',2);
grid on;
title('$\hat{J}_{yy} (kg \cdot m^2)-Estimation\,step$','Interpreter','latex')
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
title('$\hat{y}_{CG \rightarrow COM} (m)-Estimation\,step$','Interpreter','latex')
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
title('$\hat{x}_{CG \rightarrow COM} (m)-Estimation\,step$','Interpreter','latex')
xlabel('Estimation step(n)')
ylabel('$\hat{x}_{CG \rightarrow COM} (m)$','Interpreter','latex')

% J_xx_hat_vec(end)
J_yy_hat_vec(end)
% y_off_hat_vec(end)
% x_off_hat_vec(end)

inv(W_xx(2:10,:)'*W_xx(2:10,:))*W_xx(2:10,:)'*b_xx_vec(2:10)
inv(W_yy(2:end,:)'*W_yy(2:end,:))*W_yy(2:end,:)'*b_yy_vec

% A = W(49:51,:);
% b = M_x(49:51);
% 
% [V D] = eig(A'*A);
% 
% [V_obs D_obs] = CheckSVD(V,D);
% obs_index = CheckObsParam(V_obs);
% solution2 = GetSolution(V_obs,D_obs,A,b)
