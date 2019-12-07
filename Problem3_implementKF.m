%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Jake Vendl and Jack Toland
% ASEN 5044 - Statistical Estimation for Dynamical Systems
% Final Project - Orbit Determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

load('orbitdeterm_finalproj_KFdata.mat')

N = size(tvec,2);
T = max(tvec);

sensor = NaN(4,N,12);
for i = 2:N
    data = ydata{:,i};
    n = size(data,2);
    if n == 1
        station1 = data(4,1);
        sensor(:,i,station1) = data(:,1);
    elseif n == 2
        station1 = data(4,1);
        station2 = data(4,2);
        sensor(:,i,station1) = data(:,1);
        sensor(:,i,station2) = data(:,2);
    end    
end

plotsensor(tvec,sensor)


%% Implement Linearized Kalman Filter

LKF = zeros(4,N);
for k = 1:N
    
    
    
    
end

title = 'Linearized Kalman Filter State Trajectory';
filename = 'ASEN5044_FP_P3_LKF.png';
plottrajectory(tvec,LKF,title,filename);


%% Implement Extended Kalman Filter

EKF = zeros(4,N);
for k = 1:N
    
    
    
    
end

title = 'Extended Kalman Filter State Trajectory';
filename = 'ASEN5044_FP_P3_EKF.png';
plottrajectory(tvec,EKF,title,filename);


%% Implement Unscented Kalman Filter

UKF = zeros(4,N);
for k = 1:N
    
    
    
    
end

title = 'Unscented Kalman Filter State Trajectory';
filename = 'ASEN5044_FP_P3_UKF.png';
plottrajectory(tvec,UKF,title,filename);

function plotsensor(tvec,sensor,n)

    T = max(tvec);

    fig = figure('visible','on');
    set(fig,'Position',[100 100 900 600]);
    sgtitle('Sensor Data');

    subplot(3,1,1); hold on; grid on; grid minor;
    for i = 1:12
        plot(tvec,sensor(1,:,i),'*');
    end
    xlabel('time [sec]');
    ylabel('rho [km]');
    xlim([0 T]);

    subplot(3,1,2); hold on; grid on; grid minor;
    for i = 1:12
        plot(tvec,sensor(2,:,i),'*');
    end
    xlabel('time [sec]');
    ylabel('rho dot [km]');
    xlim([0 T]);

    subplot(3,1,3); hold on; grid on; grid minor;
    for i = 1:12
        plot(tvec,sensor(3,:,i),'*');
    end
    xlabel('time [sec]');
    ylabel('phi [km/s]');
    xlim([0 T]);

    saveas(fig,'ASEN5044_FP_P3_YDATA.png','png');
end

function plottrajectory(tvec,traj,title,filename)

    T = max(tvec);

    fig = figure('visible','on');
    set(fig,'Position',[100 100 900 600]);
    sgtitle(title);

    subplot(2,2,1); hold on; grid on; grid minor;
    plot(tvec,traj(1,:),'b-','LineWidth',1.5);
    xlabel('time [sec]');
    ylabel('X [km]');
    xlim([0 T]);

    subplot(2,2,2); hold on; grid on; grid minor;
    plot(tvec,traj(3,:),'b-','LineWidth',1.5);
    xlabel('time [sec]');
    ylabel('Y position [km]');
    xlim([0 T]);

    subplot(2,2,3); hold on; grid on; grid minor;
    plot(tvec,traj(2,:),'b-','LineWidth',1.5);
    xlabel('time [sec]');
    ylabel('Xdot [km/s]');
    xlim([0 T]);

    subplot(2,2,4); hold on; grid on; grid minor;
    plot(tvec,traj(4,:),'b-','LineWidth',1.5);
    xlabel('time [sec]');
    ylabel('Ydot [km/s]');
    xlim([0 T]);
    saveas(fig,filename,'png');

end