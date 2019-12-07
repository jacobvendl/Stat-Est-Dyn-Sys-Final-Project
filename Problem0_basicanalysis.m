%Jake Vendl and Jack Toland
%ASEN 5044 Final Project

clear all; close all; clc

T = 14000; % largest time value in given data

%this script handles the linearized KF

%load in data provided from Canvas
load('orbitdeterm_finalproj_KFdata.mat')

%set process noise
Q=Qtrue; clear Qtrue

%set measurement noise covariance
R=Rtrue; clear Rtrue

%set problem inputs
mu = 398600;             % km^3/s^2
r0 = 6678;               % km
rE = 6378;               % km
wE = 2*pi/86400;         % rad/s
dt = 10;                 % s
P = 2*pi*sqrt(r0^3/mu);  % s

%CODE OUTLINE: jack please check and maybe make your own? this is critical
%monte carlo to pick a dx
%simulate ground truth state using given Q
%simulate measurements using R
%run KF on the dx, saving out state estimation errors
%NEES and NIS tests

%STEP ONE  - generate input to truth model
x0 = [6678, 0, 0, r0*sqrt(mu/r0^3)]';
dx0 = [0, 0.075, 0, -0.021]'; %come back and get this with MC

%STEP TWO - simulate perturbed ground truth state using ode45
s0 = x0 + dx0;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, x_perturbed] = ode45(@(t,s)orbit_prop_func(t,s),tvec,s0,opts);
x_perturbed = x_perturbed';

[T,x_star] = ode45(@(t,s)orbit_prop_func(t,s),tvec,x0,opts);
x_star = x_star';

%STEP THREE = simulate ground truth measurements using ode45 result
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

%STEP FOUR - simulate linearized dynamics
Svq = chol(Q,'lower');
R = eye(3)*1e-6;%test R for a second
Svr = chol(R,'lower');
Gamma = [0 0; 1 0; 0 0; 0 1];
dx_lin = [];
dx=dx0;
dy_lin = zeros(36,length(T));
for t=1:length(T)
    %add to dx_lin for the given simulation run
    dx_lin = horzcat(dx_lin, dx(:,t));
    
    %propagate dx forward in time
    [F, Omega] = F_Gamma_variant(X(t),Y(t));
    dx(:,t+1) = F * dx(:,t) ;
    
    %loop through the stations and simulate linearized measurements
    for i=1:12
        if ~isnan(rho(i,t))
            H = H_variant(X(t),XD(t),Y(t),YD(t),Xs(i,t),XDs(i,t),Ys(i,t),YDs(i,t));
            dy(:,t) = H*dx(:,t);
            dy_lin(3*i-2:3*i,t) = dy(:,t);
        else
            dy_lin(3*i-2:3*i,t) = [nan nan nan]';
        end
    end
end
%add dx_lin to x_star to get simulated noisy ground truth states
x_sim = x_star + dx_lin;
y_sim = y_star + dy_lin;

x_sim_test = x_star + dx_lin;


%states versus time for nonlinear dynamics
figure; hold on;
sgtitle('States vs. Time, Full Nonlinear Dynamics Simulation')
subplot(4,1,1); hold on; grid on; grid minor;
plot(tvec,x_star(1,:),'-')
ylabel('X [km]')
subplot(4,1,2); hold on; grid on; grid minor; 
plot(tvec,x_star(2,:),'-')
ylabel('Xdot [km/s]')
subplot(4,1,3); hold on; grid on; grid minor; 
plot(tvec,x_star(3,:),'-')
ylabel('Y [km]')
subplot(4,1,4); hold on; grid on; grid minor; 
plot(tvec,x_star(4,:),'-')
ylabel('Ydot [km/s]')

%full nonlinear measurement data simulation
figure; hold on; 
sgtitle('Full Nonlinear Measurement Data Simulation')
subplot(3,1,1); hold on; grid on; grid minor; ylabel('rho^i [km]')
subplot(3,1,2); hold on; grid on; grid minor; ylabel('rhodot^i [km/s]')
subplot(3,1,3); hold on; grid on; grid minor; ylabel('\phi^i [rad]')
for i=1:12
    subplot(3,1,1); 
    plot(tvec,rho(i,:),'x')
    subplot(3,1,2);
    plot(tvec,rhoDot(i,:),'o')
    subplot(3,1,3);
    plot(tvec,phi(i,:),'o')
end

%total states vs time for linearized dynamics model
figure; hold on;
sgtitle(sprintf('States vs. Time, Linearized Approximate Dynamics Simulation \n dx =[%.4fkm %.4fkm/s %.4fkm %.4fkm/s]',dx0(1),dx0(2),dx0(3),dx0(4)))
subplot(4,1,1); hold on; grid on; grid minor;
title('x position [km]')
plot(tvec,x_sim(1,:),'-')
subplot(4,1,2); hold on; grid on; grid minor;
title('x velocity [km/s]')
plot(tvec,x_sim(2,:),'-')
subplot(4,1,3); hold on; grid on; grid minor;
title('y position [km]')
plot(tvec,x_sim(3,:),'-')
subplot(4,1,4); hold on; grid on; grid minor;
title('y velocity [km/s]')
plot(tvec,x_sim(4,:),'-')

%linearized approximate solution
figure; hold on; 
sgtitle('States vs. Time, Linearized Approximate Dynamics Simulation')
subplot(4,1,1); hold on; grid on; grid minor;
plot(tvec,dx_lin(1,:),'-')
ylabel('X [km]')
subplot(4,1,2); hold on; grid on; grid minor;
plot(tvec,dx_lin(2,:),'-')
ylabel('Xdot [km/s]')
subplot(4,1,3); hold on; grid on; grid minor;
plot(tvec,dx_lin(3,:),'-')
ylabel('Y [km]')
subplot(4,1,4); hold on; grid on; grid minor;
plot(tvec,dx_lin(4,:),'-')
ylabel('Ydot [km/s]')

%linearized approximate measurement simulation
figure; hold on;
sgtitle('Approximate Linearized Measurement Data Simulation')
subplot(3,1,1); hold on; grid on; grid minor; ylabel('rho^i [km]')
subplot(3,1,2); hold on; grid on; grid minor; ylabel('rhodot^i [km/s]')
subplot(3,1,3); hold on; grid on; grid minor; ylabel('\phi^i [rad]')
for i=1:12
    subplot(3,1,1); 
    plot(tvec,y_sim(3*i-2,:),'x')
    subplot(3,1,2);
    plot(tvec,y_sim(3*i-1,:),'o')
    subplot(3,1,3);
    plot(tvec,y_sim(3*i,:),'o')
end

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


function [F Omega] = F_Gamma_variant(X,Y)
mu = 398600;        % km^3/s^2
r0_nom = sqrt(X^2+Y^2);          % km
dt = 10;

A = [0, 1, 0, 0;
    (-mu*(r0_nom)^(-3))+(3*mu*X^2*r0_nom^(-5)), 0, 3*mu*X*Y*r0_nom^(-5), 0;
    0, 0, 0, 1;
    (3*mu*X*Y)*r0_nom^(-5), 0, (-mu*r0_nom^(-3))+(3*mu*Y^2*r0_nom^(-5)), 0];
F=eye(4) + dt*A;

Gamma = [0 0; 1 0; 0 0; 0 1];
Omega = dt*Gamma;
end
