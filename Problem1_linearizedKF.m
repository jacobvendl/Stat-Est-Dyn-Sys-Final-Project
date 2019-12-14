%Jake Vendl and Jack Toland
%ASEN 5044 Final Project

clear all; close all; clc

%stop if error and show workspace values
dbstop if error

%this script handles the linearized KF

%load in data provided from Canvas
load('orbitdeterm_finalproj_KFdata.mat')

%set problem inputs
mu = 398600;             % km^3/s^2
r0 = 6678;               % km
rE = 6378;               % km
wE = 2*pi/86400;         % rad/s
dt = 10;                 % s
P = 2*pi*sqrt(r0^3/mu);  % s

%STEP ONE  - generate input to truth model
x0 = [6678, 0, 0, r0*sqrt(mu/r0^3)]';
dx0 = [0, 0.001, 0, -0.001]';

%STEP TWO - simulate perturbed ground truth state using ode45
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, x_star] = ode45(@(t,s)orbit_prop_func(t,s),tvec,x0+dx0,opts);
x_star=x_star';

%make perturbed state for comparison
[T, x_nom] = ode45(@(t,s)orbit_prop_func(t,s),tvec,x0,opts);
x_nom=x_nom';

%STEP THREE = simulate ground truth measurements using ode45 result
X=x_nom(1,:); Y=x_nom(3,:); XD=x_nom(2,:); YD=x_nom(4,:);
Xs = zeros(12,length(T));
Ys = zeros(12,length(T));
XDs = zeros(12,length(T));
YDs = zeros(12,length(T));
rho = zeros(12,length(T));
rhoDot = zeros(12,length(T));
phi = zeros(12,length(T));
y_nom = zeros(36,length(T));
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
        y_nom(3*i-2,t) = rho(i,t);
        y_nom(3*i-1,t) = rhoDot(i,t);
        y_nom(3*i,t) = phi(i,t);
    end
end

%TMT
%set process noise
Q_process=eye(2)*1e-9;
%set measurement noise covariance
R=zeros(3); R(1,1)=0.01;R(2,2)=1;R(3,3)=0.01;
Svq = chol(Q_process,'lower');
Svr = chol(R,'lower');
Omega = [0 0; 1 0; 0 0; 0 1];


Nsim=50;
NEESamps = zeros(Nsim,length(tvec)-1);
NISamps = zeros(Nsim,length(tvec)-1);

for s=1:Nsim
    x_noisy = zeros(4,length(T));
    rho_noisy = zeros(12,length(T));
    rhoDot_noisy = zeros(12,length(T));
    phi_noisy = zeros(12,length(T));
    y_noisy = zeros(36,length(T));
    for k=1:length(tvec)
        %calculate process noise and measurement noise for the time step
        qk = randn(2,1);
        wk = (Svq*qk);
        
        %add the noise to the state output
        x_noisy(:,k) = x_star(:,k) + dt*Omega*wk;
        
        %loop through the stations and simulate linearized measurements
        for i=1:12
            rk = randn(3,1);
            vk = (Svr*rk);
            phi_noisy(i,k) = atan2((x_noisy(3,k)-Ys(i,k)),(x_noisy(1,k)-Xs(i,k)));
            thetaCheck = atan2(Ys(i,k),Xs(i,k));
            if (thetaCheck-pi/2) > (thetaCheck+pi/2)
                upperBound = thetaCheck-pi/2;
                lowerBound = thetaCheck+pi/2;
            else
                upperBound = thetaCheck+pi/2;
                lowerBound = thetaCheck-pi/2;
            end
            if (lowerBound <= phi(i,k) && phi(i,k) <= upperBound) ...
                    || (lowerBound-2*pi <= phi(i,k) && phi(i,k)<=upperBound-2*pi)... %accomodate phi wrapping
                    || (lowerBound+2*pi <= phi(i,k) && phi(i,k)<=upperBound+2*pi)
                
                rho_noisy(i,k) = sqrt((x_noisy(1,k)-Xs(i,k))^2 + (x_noisy(3,k)-Ys(i,k))^2);
                rhoDot_noisy(i,k) = ((x_noisy(1,k)-Xs(i,k))*(XD(k)-XDs(i,k)) + (x_noisy(3,k)-Ys(i,k))*(YD(k)-YDs(i,k)))...
                    / rho_noisy(i,k);
            else
                rho_noisy(i,k) = nan;
                rhoDot_noisy(i,k) = nan;
                phi_noisy(i,k)=nan;
            end
            y_noisy(3*i-2,k) = rho_noisy(i,k) + vk(1);
            y_noisy(3*i-1,k) = rhoDot_noisy(i,k) + vk(2);
            y_noisy(3*i,k) = phi_noisy(i,k) + vk(3);
        end
    end
    
    %Linearized KF
    Q_KF = eye(2)*1e-8;
    R_KF = Rtrue;
    P_plus = 1e6*eye(4);
    
    dx_hat_plus = dx0;
    dx_hat_minus = zeros(4,length(tvec));
    x_hat = zeros(4,length(tvec));
    NEESsshist = zeros(1,length(tvec)-1);
    NISsshist = zeros(1,length(tvec)-1);
    for k=1:length(tvec)-1
        
        %KF steps now that we have noisy data to work with
        [F, Omega] = F_Omega_variant(X(k),Y(k));
        dx_hat_minus(:,k+1) = F*dx_hat_plus(:,k);
        P_minus = F*P_plus*F' + Omega*Q_KF*Omega';
        
        H = [];
        dy_KF = [];
        R = R_KF;
        %loop through the stations to establish sensor measurement at k
        for i=1:12
            if ~isnan(rho(i,k+1))
                dy_KF = vertcat(dy_KF, y_noisy(3*i-2:3*i,k+1)-y_nom(3*i-2:3*i,k+1));
                H = vertcat(H,H_variant(X(k+1),XD(k+1),Y(k+1),YD(k+1),Xs(i,k+1),XDs(i,k+1),Ys(i,k+1),YDs(i,k+1)));
                if length(dy_KF) >= 4
                    R = blkdiag(R_KF,R_KF);
                end
            end
        end
        if isempty(H)==1
            K=zeros(4,3);
            H=zeros(3,4);
            dy_KF=zeros(3,1);
        else
            Sk = H*P_minus*H' + R;
            K = P_minus*H'*inv(H*P_minus*H' + R);
        end
        %Joseph formulation
        P_plus =(eye(4)-K*H)*P_minus*(eye(4)-K*H)' + K*R*K';
        
        %update state prediction
        dx_hat_plus(:,k+1) = dx_hat_minus(:,k+1) + K*(dy_KF - H*dx_hat_minus(:,k+1));
        x_hat(:,k) = x_nom(:,k) + dx_hat_plus(:,k);
        
        %save off covariance info
        twoSigX(k) = 2*sqrt(P_plus(1,1));
        twoSigXdot(k) = 2*sqrt(P_plus(2,2));
        twoSigY(k) = 2*sqrt(P_plus(3,3));
        twoSigYdot(k) = 2*sqrt(P_plus(4,4));
        
        %compute NEES and NIS statistics
        NEESsshist(k) = (x_noisy(:,k)-x_hat(:,k))'*inv(P_plus)*(x_noisy(:,k)-x_hat(:,k));
        NISsshist(k) = dy_KF'*inv(Sk)*dy_KF; %/ (length(dy_KF)/3);
    end
    NEESamps(s,:) = NEESsshist;
    NISamps(s,:) = NISsshist;
    fprintf('s=%.0f \n',s)
end
twoSigX(1401) = 2*sqrt(P_plus(1,1));
twoSigXdot(1401) = 2*sqrt(P_plus(2,2));
twoSigY(1401) = 2*sqrt(P_plus(3,3));
twoSigYdot(1401) = 2*sqrt(P_plus(4,4));

%plot NEES statistics
epsNEESbar = mean(NEESamps,1);
alphaNEES = 0.05; %significance level
Nnx = Nsim*4; %N*n
r1x = chi2inv(alphaNEES/2,Nnx)./Nsim;
r2x = chi2inv(1-alphaNEES/2,Nnx)./Nsim;
fig = figure; hold on; grid on; grid minor;
set(fig,'Position',[100 100 900 600]);
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',1.5)
plot(r1x*ones(size(epsNEESbar)),'k--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'k--','LineWidth',2)
ylabel('NEES Statistics, avg \epsilon_x')
xlabel('time step k')
title(sprintf('LKF, NEES Estimation Results, N=%.0f',Nsim))
legend('NEES @ time k','r_1 bound','r_2 bound')
ylim([0 10])
saveas(fig,'Problem1_NEES.png','png');

epsNISbar = mean(NISamps,1);
alphaNIS = 0.03; %significance level
Nny = Nsim*3; %N*p
r1y = chi2inv(alphaNIS/2,Nny)./Nsim;
r2y = chi2inv(1-alphaNIS/2,Nny)./Nsim;
fig = figure; hold on; grid on; grid minor;
set(fig,'Position',[100 100 900 600]);
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',1.5)
plot(r1y*ones(size(epsNISbar)),'k--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'k--','LineWidth',2)
ylabel('NIS Statistics, avg \epsilon_y')
xlabel('time step k')
title(sprintf('LKF, NIS Estimation Results, N=%.0f',Nsim))
legend('NIS @ time k','r_1 bound','r_2 bound')
%ylim([0 10])
saveas(fig,'Problem1_NIS.png','png');

%save out one of the simulated states for plotting purposes
x_hat = x_star + dx_hat_plus;

%total states vs time for noisy linearized dynamics model
fig = figure; hold on;
set(fig,'Position',[100 100 900 600]);
sgtitle(sprintf('LKF, Simulated Noisy Ground Truth Data \n dx =[%.4fkm %.4fkm/s %.4fkm %.4fkm/s]',dx0(1),dx0(2),dx0(3),dx0(4)))
subplot(4,1,1); hold on; grid on; grid minor;
ylabel('X [km]')
plot(tvec,x_noisy(1,:),'b-','LineWidth',1.25)
subplot(4,1,2); hold on; grid on; grid minor;
ylabel('Xdot [km/s]')
plot(tvec,x_noisy(2,:),'b-','LineWidth',1.25)
subplot(4,1,3); hold on; grid on; grid minor;
ylabel('Y [km]')
plot(tvec,x_noisy(3,:),'b-','LineWidth',1.25)
subplot(4,1,4); hold on; grid on; grid minor;
ylabel('Ydot [km/s]'); xlabel('Time [s]')
plot(tvec,x_noisy(4,:),'b-','LineWidth',1.25)
saveas(fig,'Problem1_Noisy_GT.png','png');

%noisy linearized approximate measurement simulation
fig = figure; hold on;
set(fig,'Position',[100 100 900 600]);
sgtitle('LKF, Simulated Noisy Measurement Data')
subplot(3,1,1); hold on; grid on; grid minor; ylabel('rho^i [km]')
subplot(3,1,2); hold on; grid on; grid minor; ylabel('rhodot^i [km/s]')
subplot(3,1,3); hold on; grid on; grid minor; ylabel('\phi^i [rad]')
for i=1:12
    subplot(3,1,1);
    plot(tvec,y_noisy(3*i-2,:),'x')
    subplot(3,1,2);
    plot(tvec,y_noisy(3*i-1,:),'o')
    subplot(3,1,3);
    plot(tvec,y_noisy(3*i,:),'o')
end
saveas(fig,'Problem1_Noisy_Y.png','png');



fig = figure; hold on;
set(fig,'Position',[100 100 900 600]);
sgtitle('LKF, State Estimation Errors')
subplot(4,1,1); hold on; grid on; grid minor;
plot(tvec,x_hat(1,:)-x_star(1,:),'b-','LineWidth',1.25)
plot(tvec,twoSigX,'k--','LineWidth',1)
plot(tvec,-twoSigX,'k--','LineWidth',1)
legend('x_{hat} - x_{true}','+/- 2\sigma')
ylabel('X [km]')
subplot(4,1,2); hold on; grid on; grid minor;
plot(tvec,x_hat(2,:)-x_star(2,:),'b-','LineWidth',1.25)
plot(tvec,twoSigXdot,'k--','LineWidth',1)
plot(tvec,-twoSigXdot,'k--','LineWidth',1)
ylabel('Xdot [km/s]')
ylim([-1 1])
subplot(4,1,3); hold on; grid on; grid minor;
plot(tvec,x_hat(3,:)-x_star(3,:),'b-','LineWidth',1.25)
plot(tvec,twoSigY,'k--','LineWidth',1)
plot(tvec,-twoSigY,'k--','LineWidth',1)
ylabel('Y [km]')
subplot(4,1,4); hold on; grid on; grid minor;
plot(tvec,x_hat(4,:)-x_star(4,:),'b-','LineWidth',1.25)
plot(tvec,twoSigYdot,'k--','LineWidth',1)
plot(tvec,-twoSigYdot,'k--','LineWidth',1)
ylabel('Ydot [km/s]'); xlabel('Time [s]')
ylim([-1 1])
saveas(fig,'Problem1_Error.png','png');


% figure; hold on;
% sgtitle(sprintf('Linearized KF plotted against Nonlinear Perturbed Simulation \n dx=[%.4fkm %.4fkm/s %.4fkm %.4fkm/s]',dx0(1),dx0(2),dx0(3),dx0(4)))
% subplot(4,1,1); hold on; grid on; grid minor;
% plot(tvec,x_star(1,:),'b-','LineWidth',2)
% plot(tvec,x_hat(1,:),'r--','LineWidth',2)
% legend('ode45 perturbed','xhat+')
% ylabel('X [km]')
% subplot(4,1,2); hold on; grid on; grid minor;
% plot(tvec,x_star(2,:),'b-','LineWidth',2)
% plot(tvec,x_hat(2,:),'r--','LineWidth',2)
% ylabel('Xdot [km/s]')
% subplot(4,1,3); hold on; grid on; grid minor;
% plot(tvec,x_star(3,:),'b-','LineWidth',2)
% plot(tvec,x_hat(3,:),'r--','LineWidth',2)
% ylabel('Y [km]')
% subplot(4,1,4); hold on; grid on; grid minor;
% plot(tvec,x_star(4,:),'b-','LineWidth',2)
% plot(tvec,x_hat(4,:),'r--','LineWidth',2)
% ylabel('Ydot [km/s]'); xlabel('Time [s]')



% figure; hold on
% sgtitle('Linearized Approx Perturbations vs Time')
% subplot(4,1,1); hold on; grid on; grid minor;
% plot(tvec,x_hat(1,:)-x_perturb(1,:),'b-','LineWidth',2)
% ylabel('X [km]')
% subplot(4,1,2); hold on; grid on; grid minor;
% plot(tvec,x_hat(2,:)-x_perturb(2,:),'b-','LineWidth',2)
% ylabel('Xdot [km/s]')
% subplot(4,1,3); hold on; grid on; grid minor;
% plot(tvec,x_hat(3,:)-x_perturb(3,:),'b-','LineWidth',2)
% ylabel('Y [km]')
% subplot(4,1,4); hold on; grid on; grid minor;
% plot(tvec,x_hat(4,:)-x_perturb(4,:),'b-','LineWidth',2)
% ylabel('Ydot [km/s]'); xlabel('Time [s]')



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
