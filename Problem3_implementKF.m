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

x0 = [6678, 0, 0, r0*sqrt(mu/r0^3)]';
dx0 = [0, 0.075, 0, -0.021]';

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

% Calculate sensor positions
sensor_pos = zeros(4,N,12);
for i = 1:12
    theta = (i-1)*pi/6;
    for t = 1:N
        currentTime = tvec(t);
        sensor_pos(1,t,i) = rE*cos(wE*currentTime + theta);
        sensor_pos(2,t,i) = rE*sin(wE*currentTime + theta);
        sensor_pos(3,t,i) = -rE*wE*sin(wE*currentTime + theta);
        sensor_pos(4,t,i) = rE*wE*cos(wE*currentTime + theta);
    end
end

%% Implement Linearized Kalman Filter

LKF = zeros(4,N);

% Nominal Orbit
[T, x_star] = ode45(@(t,s)orbit_prop_func(t,s),tvec,x0,opts);
x_star=x_star';

% % Nominal Measurements
% [y_star] = measurement_set2(tvec,x_star);
% 
% plotsensor(tvec,y_star)

% Initialized KF
P_plus = 1e3*eye(4);
dx_hat_plus = dx0;
dx_hat_minus = zeros(4,N);
x_hat = zeros(4,N);

Q_KF=eye(2)*1e-9;

for k = 1:N-1
    [F, Gamma] = F_Gamma_variant(x_star(1,k),x_star(3,k));
    dx_hat_minus(:,k+1) = F*dx_hat_plus(:,k);
    P_minus = F*P_plus*F' + Gamma*Q_KF*Gamma';
    
    H = [];
    dy_KF = [];
    R_KF = R;
    for i = 1:12
        if ~isnan(sensor(1,k,i))
            y_star = measurement(x_star(:,k),sensor_pos(:,k+1,i));
            dy_KF = vertcat(dy_KF, sensor(1:3,k,i)-y_star);
            H = vertcat(H,H_variant(x_star(:,k+1),sensor_pos(:,k+1,i)));
            if length(dy_KF) >= 4
                R_KF = blkdiag(R,R);
            end
        end
    end
    if isempty(H) == 1
        K = zeros(4,3);
        H = zeros(3,4);
        dy_KF = zeros(3,1);
    else
        Sk = H*P_minus*H' + R_KF;
        K = P_minus*H'*inv(H*P_minus*H' + R_KF);
    end
    P_plus = (eye(4) - K*H)*P_minus;
    
    dx_hat_plus(:,k+1) = dx_hat_minus(:,k+1) + K*(dy_KF - H*dx_hat_minus(:,k+1));
    x_hat(:,k) = x_star(:,k) + dx_hat_plus(:,k);
    LKF(:,k) = x_hat(:,k);
   
    %save off covariance info
    twoSigX(k) = 2*sqrt(P_plus(1,1));
    twoSigXdot(k) = 2*sqrt(P_plus(2,2));
    twoSigY(k) = 2*sqrt(P_plus(3,3));
    twoSigYdot(k) = 2*sqrt(P_plus(4,4));
end

title = 'Linearized Kalman Filter State Trajectory';
filename = 'ASEN5044_FP_P3_LKF.png';
plottrajectory(tvec,LKF,title,filename);


%% Implement Extended Kalman Filter

% Step 1: Initialize with x_plus_0 and P_plus_0


EKF = zeros(4,N);
x_plus = x0;
P_plus = eye(4)*1e-6;


for k = 2:N-1 % k represents k+1
    
    ts = tvec(k-1);
    tf = tvec(k);
   
    [F, Gamma] = F_Gamma_variant(x_plus(1),x_plus(3));
    
    
    % 1) Time update for k+1
    [~, x] = ode45(@(t,s)orbit_prop_func(t,s),[ts tf],x_plus,opts);
    x_minus = x(end,:)';
    P_minus = F*P_plus*F' + Gamma*Q*Gamma';
    
    
    % 2) Measurement Update for k+1, work through all stations
    H = [];
    e_KF = [];
    R_KF = [];
    for i = 1:12
        if ~isnan(sensor(1,k,i))
            ynom_minus = measurement(x_minus,sensor_pos(:,k,i));
            H = H_variant(x_minus,sensor_pos(:,k,i));
            
            e = sensor(1:3,k,i) - ynom_minus;
            e_KF = vertcat(e_KF,e);
            
            if length(e_KF) >= 4
                R_KF = blkdiag(R,R);
            end
            
            
        end
    end
    if isempty(H) == 1
        K = zeros(4,3);
        H = zeros(3,4);
        e_KF = zeros(3,1);
    else
        K = P_minus*H'*inv(H*P_minus*H'+R);
    end    
    x_plus = x_minus + K*e;
    P_plus = (eye(4) - K*H)*P_minus;
    
    
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

function [y_star] = measurement_set2(T,x_star)

rE = 6378;               % km
wE = 2*pi/86400;         % rad/s

X=x_star(1,:); Y=x_star(3,:); XD=x_star(2,:); YD=x_star(4,:);
Xs = zeros(12,length(T));
Ys = zeros(12,length(T));
XDs = zeros(12,length(T));
YDs = zeros(12,length(T));
rho = zeros(12,length(T));
rhoDot = zeros(12,length(T));
phi = zeros(12,length(T));
num = zeros(12,length(T));
y_star = zeros(3,length(T),12);
%now simulate the measurements for all time
for i=1:12 %stations
    theta = (i-1)*pi/6;
    for t=1:length(T) %loop through one orbit period
        currentTime = T(t);
        
        %find station position and velocity
        Xs(i,t) = rE*cos(wE*currentTime + theta);
        Ys(i,t) = rE*sin(wE*currentTime + theta);
        XDs(i,t) = -rE*wE*sin(wE*currentTime + theta);
        YDs(i,t) = rE*wE*cos(wE*currentTime + theta);
        
        %peform check at given time to see if s/c is visible
        phi(i,t) = atan2((Y(t)-Ys(i,t)),(X(t)-Xs(i,t)));
        thetaCheck = atan2(Ys(i,t),Xs(i,t));
        if (thetaCheck-pi/2) > (thetaCheck+pi/2)
            upperBound = thetaCheck-pi/2;
            lowerBound = thetaCheck+pi/2;
        else
            upperBound = thetaCheck+pi/2;
            lowerBound = thetaCheck-pi/2;
        end
        if (lowerBound <= phi(i,t) && phi(i,t) <= upperBound) ...
                || (lowerBound-2*pi <= phi(i,t) && phi(i,t)<=upperBound-2*pi)... %accomodate phi wrapping
                || (lowerBound+2*pi <= phi(i,t) && phi(i,t)<=upperBound+2*pi)
            
            rho(i,t) = sqrt((X(t)-Xs(i,t))^2 + (Y(t)-Ys(i,t))^2);
            rhoDot(i,t) = ((X(t)-Xs(i,t))*(XD(t)-XDs(i,t)) + (Y(t)-Ys(i,t))*(YD(t)-YDs(i,t)))...
                / rho(i,t);
        else
            rho(i,t) = nan;
            rhoDot(i,t) = nan;
            phi(i,t)=nan;
            num(i,t)=nan;
        end
        y_star(1,:,i) = rho(i,t);
        y_star(2,:,i) = rhoDot(i,t);
        y_star(3,:,i) = phi(i,t);
        y_star(4,:,i) = num(i,t);
    end
end
end


function [y] = measurement(state,station)  
y = zeros(3,1);

% rho
y(1) = sqrt((state(1)-station(1))^2 + (state(3)-station(3))^2);
% rhodot
y(2) = ((state(1)-station(1))*(state(2)-station(2))...
    + (state(3)-station(3))*(state(4)-station(4)))/ y(1);
% phi
y(3) = atan2((state(3)-station(3)),(state(1)-station(1)));
end


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

function [ H ] = H_variant(state, station)
X = state(1);
Y = state(2);
Xdot = state(3);
Ydot = state(4);

Xs = station(1);
Ys = station(2);
Xsdot = station(3);
Ysdot = station(4);

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


function [F Gamma] = F_Gamma_variant(X,Y)
mu = 398600;        % km^3/s^2
r0_nom = 6678;          % km
dt = 10;

A = [0, 1, 0, 0;
    (-mu*(r0_nom)^(-3))+(3*mu*X^2*r0_nom^(-5)), 0, 3*mu*X*Y*r0_nom^(-5), 0;
    0, 0, 0, 1;
    (3*mu*X*Y)*r0_nom^(-5), 0, (-mu*r0_nom^(-3))+(3*mu*Y^2*r0_nom^(-5)), 0];
F=eye(4) + dt*A;

Gamma = [0 0; 1 0; 0 0; 0 1];
Gamma = dt*Gamma;
end

