%% motor calibration
%loads and plots the motor calibration data
%path and file name of data
figure
fpath = 'C:\Users\etuthill\OneDrive - Olin College of Engineering\ESA\Rocky_Project\rocky_project\'; %path (change this!)
fname_in = 'motor_data.txt'; %file name (change this!)
%load the motor calibration data
motor_data = importdata([fpath,fname_in]);
%unpack the motor calibration data
t_m = motor_data(830:1800,1);
y_L = motor_data(830:1800,2); v_L = motor_data(830:1800,3);
y_R = motor_data(830:1800,4); v_R = motor_data(830:1800,5);
%plot the motor calibration data
figure(1);
subplot(2,1,1);
hold on
plot(t_m,v_L,'k','linewidth',1);
plot(t_m,v_R,'r','linewidth',1);
xlabel('time (sec)'); ylabel('wheel speed (m/sec)');
title('Motor Calibration Data');
xlim([15, 26])
h1 = legend('Left Wheel','Right Wheel');
set(h1,'location','bestoutside');
subplot(2,1,2);
hold on
plot(t_m,y_L,'k','linewidth',1);
plot(t_m,y_R,'r--','linewidth',1);
xlim([15, 26])
xlabel('time (sec)'); ylabel('wheel command (-)');
title('Motor Calibration Data');
h2 = legend('Left Wheel','Right Wheel');
set(h2,'location','bestoutside');
hold off



%% determine alpha and beta

t_fit = t_m(21:120);
v_fit = v_L(21:120);
y_fit = y_L(21:120);

t_fit = t_fit - t_fit(1);
ft = fittype('c*(1-exp(-a*x))');
model = fit(t_fit,v_fit,ft);

t_fit - t_fit - t_fit(1);

a = model.a;
c = model.c;

Y = mean(y_fit(end-10:end));

alpha = a;
beta = c/Y;

alpha
beta

%% gyro calibration
figure
%path and file name of data
fpath = 'C:\Users\etuthill\OneDrive - Olin College of Engineering\ESA\Rocky_Project\rocky_project\'; %path (change this!)
fname_in = 'gyro_data.txt'; %file name (change this!)
%load the pendulum calibration data
pendulum_data = importdata([fpath,fname_in]);
%unpack the pendulum calibration data
t_g = pendulum_data(1055:end,1); theta = pendulum_data(1055:end,2);
%plot the motor calibration data
figure(1);
hold on
plot(t_g,theta,'k','linewidth',1);
xlabel('time (sec)'); ylabel('angle (rad)');
title('Pendulum Calibration Data');

%% natural frequency

t_fit = t_g - t_g(1);

t_fit = t_fit(:);
x = theta(:);

signal_params = regress_underdamped_response(t_fit,x);

omega_n = signal_params.omega_n;
zeta = signal_params.zeta;
f_n = omega_n/(2*pi);

omega_n_2 = omega_n^2

%% effective length

l_eff = 9.81/omega_n_2