%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jake Vendl and Jack Toland
% ASEN 5044 - Statistical Estimation for Dynamical Systems
% Final Project - Orbit Determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%load in provided data
load('orbitdeterm_finalproj_KFdata.mat')

%set problem inputs
mu = 398600;             % km^3/s^2
r0 = 6678;               % km
rE = 6378;               % km
wE = 2*pi/86400;         % rad/s
dt = 10;                 % s
P = 2*pi*sqrt(r0^3/mu);  % s

x0 = [6678, 0, 0, r0*sqrt(mu/r0^3)]';
dx0 = [0, 0.01, 0, 0.01]';

%generate truth reference trajectory, unperturbed
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, x_star] = ode45(@(t,s)orbit_prop_func(t,s),tvec,x0,opts);
x_star=x_star';

%generate truth data y_star
X=x_star(1,:); Y=x_star(3,:); XD=x_star(2,:); YD=x_star(4,:);
Xs = zeros(12,length(T));
Ys = zeros(12,length(T));
XDs = zeros(12,length(T));
YDs = zeros(12,length(T));
rho = zeros(12,length(T));
rhoDot = zeros(12,length(T));
phi = zeros(12,length(T));
y_star = zeros(36,length(T));
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
        end
        y_star(3*i-2,t) = rho(i,t);
        y_star(3*i-1,t) = rhoDot(i,t);
        y_star(3*i,t) = phi(i,t);
    end
end

%reformat incoming data to reflect the format of y_star
y_noisy = zeros(36,length(tvec));
for t=1:length(tvec)
    for i=1:12
        s = size(ydata{t});
        if s(2) == 1 %there was one measurement recorded at the time
            if ydata{t}(4,1) == i
                y_noisy(3*i-2:3*i,t) = ydata{t}(1:3);
            else
                y_noisy(3*i-2:3*i,t) = [NaN;NaN;NaN]';
            end
        elseif s(2) == 2 %there were two measurements recorded
            if ydata{t}(4,1) == i
                y_noisy(3*i-2:3*i,t) = ydata{t}(1:3);
            elseif ydata{t}(4,2) == i
                y_noisy(3*i-2:3*i,t) = ydata{t}(1:3);
            else
                y_noisy(3*i-2:3*i,t) = [NaN;NaN;NaN]';
            end
        else
            y_noisy(3*i-2:3*i,t) = [NaN;NaN;NaN]';
        end
    end
end


Q_KF = Qtrue;
R = Rtrue;
Q_KF = eye(2)*1e-10;

% Implement Linearized Kalman Filter
LKF = zeros(4,length(tvec));

P_plus = eye(4)*1e4;
dx_hat_plus = dx0;
dx_hat_minus = zeros(4,length(tvec));
for k=1:length(tvec)-1
    [F, Omega] = F_Omega_variant(X(k),Y(k));
    dx_hat_minus(:,k+1) = F*dx_hat_plus(:,k);
    P_minus = F*P_plus*F' + Omega*Q_KF*Omega';
    
    H = [];
    innov = [];
    R_KF = R;
    %find innov based on provided data
    for i=1:12
        if ~isnan(y_noisy(3*i,k+1)) && ~isnan(y_star(3*i,k+1))
            innov = vertcat(innov, y_noisy(3*i-2:3*i,k+1)-y_star(3*i-2:3*i,k+1));
            H = vertcat(H,H_variant(X(k+1),XD(k+1),Y(k+1),YD(k+1),Xs(i,k+1),XDs(i,k+1),Ys(i,k+1),YDs(i,k+1)));
            if length(innov) >= 4
                R_KF = blkdiag(R,R);
            end
        end
    end
    if isempty(H)==1
        K=zeros(4,3);
        H=zeros(3,4);
        innov=zeros(3,1);
    else
        Sk = H*P_minus*H' + R_KF;
        K = P_minus*H'*inv(H*P_minus*H' + R_KF);
    end
    P_plus = (eye(4) - K*H)*P_minus;
    
    %update state prediction
    dx_hat_plus(:,k+1) = dx_hat_minus(:,k+1) + K*(innov - H*dx_hat_minus(:,k+1));
    x_hat(:,k) = x_star(:,k) + dx_hat_plus(:,k);
    
    %save off covariance info
    twoSigX(k) = 2*sqrt(P_plus(1,1));
    twoSigXdot(k) = 2*sqrt(P_plus(2,2));
    twoSigY(k) = 2*sqrt(P_plus(3,3));
    twoSigYdot(k) = 2*sqrt(P_plus(4,4));
end
x_hat(:,1401) = x_star(:,1401) + dx_hat_plus(:,1401);
twoSigX(1401) = 2*sqrt(P_plus(1,1));
twoSigXdot(1401) = 2*sqrt(P_plus(2,2));
twoSigY(1401) = 2*sqrt(P_plus(3,3));
twoSigYdot(1401) = 2*sqrt(P_plus(4,4));

figure; hold on; grid on; grid minor;
sgtitle('LKF Estimated States for provided data')
subplot(4,1,1); hold on; grid on; grid minor
plot(tvec,x_hat(1,:),'b-','LineWidth',2)
plot(tvec,x_hat(1,:) + twoSigX,'k--','LineWidth',1)
plot(tvec,x_hat(1,:) - twoSigX,'k--','LineWidth',1)
ylabel('X [km]')
legend('State','+/- 2\sigma')
subplot(4,1,2); hold on; grid on; grid minor
plot(tvec,x_hat(2,:),'b-','LineWidth',2)
plot(tvec,x_hat(2,:) + twoSigXdot,'k--','LineWidth',1)
plot(tvec,x_hat(2,:) - twoSigXdot,'k--','LineWidth',1)
ylabel('XDot [km/s]')
subplot(4,1,3); hold on; grid on; grid minor
plot(tvec,x_hat(3,:),'b-','LineWidth',2)
plot(tvec,x_hat(3,:) + twoSigY,'k--','LineWidth',1)
plot(tvec,x_hat(3,:) - twoSigY,'k--','LineWidth',1)
ylabel('Y [km]')
subplot(4,1,4); hold on; grid on; grid minor
plot(tvec,x_hat(4,:),'b-','LineWidth',2)
plot(tvec,x_hat(4,:) + twoSigYdot,'k--','LineWidth',1)
plot(tvec,x_hat(4,:) - twoSigYdot,'k--','LineWidth',1)
xlabel('time [s]')
ylabel('YDot [km/s]')

figure; hold on;
sgtitle('LKF - State Estimation Errors')
subplot(4,1,1); hold on; grid on; grid minor;
plot(tvec,x_hat(1,:)-x_star(1,:),'b-','LineWidth',2)
legend('x_{hat} - x*')
ylabel('X [km]')
subplot(4,1,2); hold on; grid on; grid minor;
plot(tvec,x_hat(2,:)-x_star(2,:),'b-','LineWidth',2)
ylabel('Xdot [km/s]')
subplot(4,1,3); hold on; grid on; grid minor;
plot(tvec,x_hat(3,:)-x_star(3,:),'b-','LineWidth',2)
ylabel('Y [km]')
subplot(4,1,4); hold on; grid on; grid minor;
plot(tvec,x_hat(4,:)-x_star(4,:),'b-','LineWidth',2)
ylabel('Ydot [km/s]')




%propagation function
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


function [F Omega] = F_Omega_variant(X,Y)
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
