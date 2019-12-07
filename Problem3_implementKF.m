%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Jake Vendl and Jack Toland
% ASEN 5044 - Statistical Estimation for Dynamical Systems
% Final Project - Orbit Determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

load('orbitdeterm_finalproj_KFdata.mat')

N = size(tvec,2);


%% Implement Linearized Kalman Filter

LKF = zeros(4,N);
for k = 1:N
    
    
    
    
end

fig = figure('visible','on');
set(fig,'Position',[100 100 900 600]);
sgtitle('Extended Kalman Filter State Trajectory');

subplot(2,2,1); hold on; grid on; grid minor;
plot(tvec,LKF(:,1),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('X [km]');
xlim([0 T]);

subplot(2,2,2); hold on; grid on; grid minor;
plot(tvec,LKF(:,3),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Y position [km]');
xlim([0 T]);

subplot(2,2,3); hold on; grid on; grid minor;
plot(tvec,LKF(:,2),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Xdot [km/s]');
xlim([0 T]);

subplot(2,2,4); hold on; grid on; grid minor;
plot(tvec,LKF(:,4),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Ydot [km/s]');
xlim([0 T]);
saveas(fig,'ASEN5044_FP_P3_LKF.png','png');


%% Implement Extended Kalman Filter

EKF = zeros(4,N);
for k = 1:N
    
    
    
    
end

fig = figure('visible','on');
set(fig,'Position',[100 100 900 600]);
sgtitle('Extended Kalman Filter State Trajectory');

subplot(2,2,1); hold on; grid on; grid minor;
plot(tvec,EKF(:,1),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('X [km]');
xlim([0 T]);

subplot(2,2,2); hold on; grid on; grid minor;
plot(tvec,EKF(:,3),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Y position [km]');
xlim([0 T]);

subplot(2,2,3); hold on; grid on; grid minor;
plot(tvec,EKF(:,2),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Xdot [km/s]');
xlim([0 T]);

subplot(2,2,4); hold on; grid on; grid minor;
plot(tvec,EKF(:,4),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Ydot [km/s]');
xlim([0 T]);
saveas(fig,'ASEN5044_FP_P3_EKF.png','png');


%% Implement Unscented Kalman Filter

UKF = zeros(4,N);
for k = 1:N
    
    
    
    
end

fig = figure('visible','on');
set(fig,'Position',[100 100 900 600]);
sgtitle('Unscented Kalman Filter State Trajectory');

subplot(2,2,1); hold on; grid on; grid minor;
plot(tvec,UKF(:,1),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('X [km]');
xlim([0 T]);

subplot(2,2,2); hold on; grid on; grid minor;
plot(tvec,UKF(:,3),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Y position [km]');
xlim([0 T]);

subplot(2,2,3); hold on; grid on; grid minor;
plot(tvec,UKF(:,2),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Xdot [km/s]');
xlim([0 T]);

subplot(2,2,4); hold on; grid on; grid minor;
plot(tvec,UKF(:,4),'b-','LineWidth',1.5);
xlabel('time [sec]');
ylabel('Ydot [km/s]');
xlim([0 T]);
saveas(fig,'ASEN5044_FP_P3_UKF.png','png');