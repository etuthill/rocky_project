%% motor calibration
%loads and plots the motor calibration data
%path and file name of data
figure
fpath = 'C:\Users\etuthill\OneDrive - Olin College of Engineering\ESA\Rocky_Project\'; %path (change this!)
fname_in = 'motor_data.txt'; %file name (change this!)
%load the motor calibration data
motor_data = importdata([fpath,fname_in]);
%unpack the motor calibration data
t = motor_data(830:1800,1);
y_L = motor_data(830:1800,2); v_L = motor_data(830:1800,3);
y_R = motor_data(830:1800,4); v_R = motor_data(830:1800,5);
%plot the motor calibration data
figure(1);
subplot(2,1,1);
hold on
plot(t,v_L,'k','linewidth',1);
plot(t,v_R,'r','linewidth',1);
xlabel('time (sec)'); ylabel('wheel speed (m/sec)');
title('Motor Calibration Data');
xlim([15, 26])
h1 = legend('Left Wheel','Right Wheel');
set(h1,'location','bestoutside');
subplot(2,1,2);
hold on
plot(t,y_L,'k','linewidth',1);
plot(t,y_R,'r--','linewidth',1);
xlim([15, 26])
xlabel('time (sec)'); ylabel('wheel command (-)');
title('Motor Calibration Data');
h2 = legend('Left Wheel','Right Wheel');
set(h2,'location','bestoutside');
hold off

%% gyro calibration
figure
%path and file name of data
fpath = 'C:\Users\etuthill\OneDrive - Olin College of Engineering\ESA\Rocky_Project\'; %path (change this!)
fname_in = 'gyro_data.txt'; %file name (change this!)
%load the pendulum calibration data
pendulum_data = importdata([fpath,fname_in]);
%unpack the pendulum calibration data
t = pendulum_data(1055:end,1); theta = pendulum_data(1055:end,2);
%plot the motor calibration data
figure(1);
hold on
plot(t,theta,'k','linewidth',1);
xlabel('time (sec)'); ylabel('angle (rad)');
title('Pendulum Calibration Data');