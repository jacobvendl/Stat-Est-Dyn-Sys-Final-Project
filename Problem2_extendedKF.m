%Jake Vendl and Jack Toland

close all; clear all; clc

%set problem inputs
mu = 398600;             % km^3/s^2
r0 = 6678;               % km
rE = 6378;               % km
wE = 2*pi/86400;         % rad/s
dt = 10;                 % s
P = 2*pi*sqrt(r0^3/mu);  % s

%load in data provided from Canvas
load('orbitdeterm_finalproj_KFdata.mat')

x0 = [6678, 0, 0, r0*sqrt(mu/r0^3)]';
dx0 = [0, 0.075, 0, -0.021]';
Gamma = [0 0; 1 0; 0 0 ; 0 1];

%set process noise
Q_KF=eye(2)*1e-16;
%set measurement noise covariance
R=zeros(3); R(1,1)=0.01;R(2,2)=1;R(3,3)=0.01;
Svq = chol(Q_KF,'lower');
Svr = chol(R,'lower');

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, x_star] = ode45(@(t,s)orbit_prop_func(t,s),tvec,x0,opts);
x_star=x_star';

%generate noisy simulated measurements
for k=1:length(tvec)
    %calculate process noise and measurement noise for the time step
    qk = randn(2,1);
    wk = (Svq*qk);
    
    %add the noise to the state output
    x_noisy = x_star(:,k) + dt*Gamma*wk;
    
    %loop through the stations and simulate linearized measurements
    for i=1:12
        theta = (i-1)*pi/6;
        currentTime = tvec(k);
        
        %find station position and velocity
        Xs = rE*cos(wE*currentTime + theta);
        Ys = rE*sin(wE*currentTime + theta);
        XDs = -rE*wE*sin(wE*currentTime + theta);
        YDs = rE*wE*cos(wE*currentTime + theta);
        
        rk = randn(3,1);
        vk = (Svr*rk);
        phi_noisy = atan2((x_noisy(3)-Ys),(x_noisy(1)-Xs));
        thetaCheck = atan2(Ys,Xs);
        if (thetaCheck-pi/2) > (thetaCheck+pi/2)
            upperBound = thetaCheck-pi/2;
            lowerBound = thetaCheck+pi/2;
        else
            upperBound = thetaCheck+pi/2;
            lowerBound = thetaCheck-pi/2;
        end
        if (lowerBound <= phi_noisy && phi_noisy <= upperBound) ...
                || (lowerBound-2*pi <= phi_noisy && phi_noisy<=upperBound-2*pi)... %accomodate phi wrapping
                || (lowerBound+2*pi <= phi_noisy && phi_noisy<=upperBound+2*pi)
            
            rho_noisy = sqrt((x_noisy(1)-Xs)^2 + (x_noisy(3)-Ys)^2);
            rhoDot_noisy = ((x_noisy(1)-Xs)*(x_noisy(2)-XDs) + (x_noisy(3)-Ys)*(x_noisy(4)-YDs)) / rho_noisy;
        else
            rho_noisy = nan;
            rhoDot_noisy = nan;
            phi_noisy=nan;
        end
        y_noisy(3*i-2:3*i,k) = [rho_noisy;rhoDot_noisy;phi_noisy] + vk;
    end
end



%Extended KF
x_hat_plus(:,1) = x0+dx0;
P_plus = eye(4)*1e4;
Q = 1e-9*eye(2);
for k=1:length(tvec)
    %call ode45 to propagate from k to k+1
    ts = tvec(k);
    tf = tvec(k+1);
    [~,temp] = ode45(@(t,s)orbit_prop_func(t,s),[ts tf],x_hat_plus(:,k),opts);
    x_hat_minus(:,k+1) = temp(end,:);
    [F,Omega] = F_Gamma_variant(x_hat_plus(1,k),x_hat_plus(3,k));
    P_minus = F*P_plus*F' + Omega*Q*Omega';
    
    %now find measurement y_hat_minus at k+1 using x_hat_minus at k+1
    X=x_hat_minus(1,k+1);
    Y=x_hat_minus(3,k+1);
    XD=x_hat_minus(2,k+1);
    YD=x_hat_minus(4,k+1);
    for i=1:12
        theta = (i-1)*pi/6;
        currentTime = tvec(k+1);
        Xs = rE*cos(wE*currentTime + theta);
        Ys = rE*sin(wE*currentTime + theta);
        XDs = -rE*wE*sin(wE*currentTime + theta);
        YDs = rE*wE*cos(wE*currentTime + theta);
        phi = atan2((Y-Ys),(X-Xs));
        thetaCheck = atan2(Ys,Xs);
        if (thetaCheck-pi/2) > (thetaCheck+pi/2)
            upperBound = thetaCheck-pi/2;
            lowerBound = thetaCheck+pi/2;
        else
            upperBound = thetaCheck+pi/2;
            lowerBound = thetaCheck-pi/2;
        end
        if (lowerBound <= phi && phi <= upperBound) ...
                || (lowerBound-2*pi <= phi && phi<=upperBound-2*pi)... %accomodate phi wrapping
                || (lowerBound+2*pi <= phi && phi<=upperBound+2*pi)
            
            rho= sqrt((X-Xs)^2 + (Y-Ys)^2);
            rhoDot = ((X-Xs)*(XD-XDs) + (Y-Ys)*(YD-YDs)) / rho;
        else
            rho = nan;
            rhoDot = nan;
            phi=nan;
        end
        y_hat_minus(:,k+1) = [rho;rhoDot;phi];
    end
    
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



