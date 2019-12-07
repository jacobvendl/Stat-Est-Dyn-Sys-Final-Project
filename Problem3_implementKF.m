%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Jake Vendl and Jack Toland
% ASEN 5044 - Statistical Estimation for Dynamical Systems
% Final Project - Orbit Determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%set problem inputs
mu = 398600;             % km^3/s^2
r0 = 6678;               % km
rE = 6378;               % km
wE = 2*pi/86400;         % rad/s
dt = 10;                 % s
P = 2*pi*sqrt(r0^3/mu);  % s

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

load('orbitdeterm_finalproj_KFdata.mat')

Q = Qtrue;
R = Rtrue;

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


x0 = [6678, 0, 0, r0*sqrt(mu/r0^3)]';

%% Implement Linearized Kalman Filter

LKF = zeros(4,N);
for k = 1:N
    
    
    
    
end

title = 'Linearized Kalman Filter State Trajectory';
filename = 'ASEN5044_FP_P3_LKF.png';
plottrajectory(tvec,LKF,title,filename);


%% Implement Extended Kalman Filter

% Step 1: Initialize with x_plus_0 and P_plus_0


EKF = zeros(4,N);
x_plus = x0;
P_plus = eye(4)*1e-6;

[F Gamma] = F_Gamma_variant(x_plus(1),x_plus(3));


for k = 1:N % k represents k+1
   
    % 1) Time update for k+1
    [~, x] = ode45(@(t,s)orbit_prop_func(t,s),[tvec(k) tvec(k+1)],x_plus,opts);
    x_minus = x(end,:)';
    P_minus = F*P_plus*F' + Gamma*Q*Gamma';
    
    
    % 2) Measurement Update for k+1
    H = H_variant(x_minus);
    ynom_minus = H*x_minus;
    
    e = sensor(:,k,1) - ynom_minus;
    K = P_minus*H'*inv(H*P_minu*H'+R);
    
    x_plus = x_minus + K*e;
    P_plus = (eye(4) - K*H)*P_minus;
    
    [F Gamma] = F_Gamma_variant(x_plus(1),x_plus(3));
    
    EKF(:,k) = x_plus;
    
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

function [ ds ] = orbit_prop_func(t,s)

mu = 398600;

x = s(1);
y = s(3);
r = sqrt(x^2+y^2);

xdot = s(2);
ydot = s(4);

xddot = -mu/r^3 * x;
yddot = -mu/r^3 * y;

ds = [xdot, xddot, ydot, yddot]';
end

function [ H ] = H_variant(X,Xdot,Y,Ydot,Xs,Xsdot,Ys,Ysdot)
%initialize H
H = zeros(3,4);

%first row
H(1,1) = (X-Xs)/sqrt((X-Xs)^2+(Y-Ys)^2);
H(1,2) = 0;
H(1,3) = (Y-Ys)/sqrt((X-Xs)^2+(Y-Ys)^2);
H(1,4) = 0;

%second row
H(2,1) = (Xdot-Xsdot)/sqrt((X-Xs)^2+(Y-Ys)^2) - (X-Xs)*((X-Xs)*(Xdot-Xsdot)+(Y-Ys)*(Ydot-Ysdot)) / ((X-Xs)^2+(Y-Ys)^2)^(3/2);
H(2,3) = (Ydot-Ysdot)/sqrt((X-Xs)^2+(Y-Ys)^2) - (Y-Ys)*((X-Xs)*(Xdot-Xsdot)+(Y-Ys)*(Ydot-Ysdot)) / ((X-Xs)^2+(Y-Ys)^2)^(3/2);
H(2,2) = (X-Xs)/sqrt((X-Xs)^2+(Y-Ys)^2);
H(2,4) = (Y-Ys)/sqrt((X-Xs)^2+(Y-Ys)^2);

%third row
H(3,1) = ((-Y+Ys)/((X-Xs)^2+(Y-Ys)^2));
H(3,2) = 0;
H(3,3) = ((X-Xs)/((X-Xs)^2+(Y-Ys)^2));
H(3,4) = 0;
end


function [F Omega] = F_Gamma_variant(X,Y)
mu = 398600;        % km^3/s^2
r0_nom = 6678;          % km
dt = 10;

A = [0, 1, 0, 0;
    (-mu*(r0_nom)^(-3))+(3*mu*X^2*r0_nom^(-5)), 0, 3*mu*X*Y*r0_nom^(-5), 0;
    0, 0, 0, 1;
    (3*mu*X*Y)*r0_nom^(-5), 0, (-mu*r0_nom^(-3))+(3*mu*Y^2*r0_nom^(-5)), 0];
F=eye(4) + dt*A;

Gamma = [0 0; 1 0; 0 0; 0 1];
Omega = dt*Gamma;
end