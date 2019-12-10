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
dx0 = [0, 0.01, 0, 0.01]';
Gamma = [0 0; 1 0; 0 0 ; 0 1];

%set process noise
Q_KF=eye(2)*1e-7; %Q_KF(1,2)=1e-12; Q_KF(2,1)=1e-12;
%set measurement noise covariance
R=zeros(3); R(1,1)=0.01;R(2,2)=1;R(3,3)=0.01;
R=eye(3)*1e-3;
Svq = chol(Q_KF,'lower');
Svr = chol(R,'lower');

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, x_star] = ode45(@(t,s)orbit_prop_func(t,s),tvec,x0,opts);
x_star=x_star';

Nsim=50;
NEESamps = zeros(Nsim,length(tvec)-1);
NISamps = zeros(Nsim,length(tvec)-1);

for s=1:Nsim
    %generate noisy simulated measurements
    for k=1:length(tvec)
        %calculate process noise and measurement noise for the time step
        qk = randn(2,1);
        wk = (Svq*qk);
        
        %add the noise to the state output
        x_noisy(:,k) = x_star(:,k) + dt*Gamma*wk;
        
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
            phi_noisy = atan2((x_noisy(3,k)-Ys),(x_noisy(1,k)-Xs));
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
                
                rho_noisy = sqrt((x_noisy(1,k)-Xs)^2 + (x_noisy(3,k)-Ys)^2);
                rhoDot_noisy = ((x_noisy(1,k)-Xs)*(x_noisy(2,k)-XDs) + (x_noisy(3,k)-Ys)*(x_noisy(4,k)-YDs)) / rho_noisy;
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
    P_plus = eye(4)*1e9;
    for k=1:length(tvec)-1
        %call ode45 to propagate from k to k+1
        ts = tvec(k);
        tf = tvec(k+1);
        [~,temp] = ode45(@(t,s)orbit_prop_func(t,s),[ts tf],x_hat_plus(:,k),opts);
        x_hat_minus(:,k+1) = temp(end,:);
        [F,Omega] = F_Omega_variant(x_hat_plus(1,k),x_hat_plus(3,k));
        P_minus = F*P_plus*F' + Omega*Q_KF*Omega';
        
        %now find measurement y_hat_minus at k+1 using x_hat_minus at k+1
        X=x_hat_minus(1,k+1);
        Y=x_hat_minus(3,k+1);
        XD=x_hat_minus(2,k+1);
        YD=x_hat_minus(4,k+1);
        y_hat_minus = [];
        H = [];
        R_KF = R;
        for i=1:12
            if ~isnan(y_noisy(3*i,k+1))
                theta = (i-1)*pi/6;
                currentTime = tvec(k+1);
                Xs = rE*cos(wE*currentTime + theta);
                Ys = rE*sin(wE*currentTime + theta);
                XDs = -rE*wE*sin(wE*currentTime + theta);
                YDs = rE*wE*cos(wE*currentTime + theta);
                phi = atan2((Y-Ys),(X-Xs));
                
                rho= sqrt((X-Xs)^2 + (Y-Ys)^2);
                rhoDot = ((X-Xs)*(XD-XDs) + (Y-Ys)*(YD-YDs)) / rho;
                y_hat_minus = vertcat(y_hat_minus,[rho;rhoDot;phi]);
                H = vertcat(H,H_variant(X,XD,Y,YD,Xs,XDs,Ys,YDs));
                if length(y_hat_minus) >= 4
                    R_KF = blkdiag(R,R);
                end
            end  
        end
        if isempty(H)==1
            K = zeros(4,3);
            H = zeros(3,4);
            innov = zeros(3,1);
        else
            %pull noisy measurement out of y_noisy matrix
            y_actual = y_noisy(~isnan(y_noisy(:,k+1)),k+1);
            innov = y_actual - y_hat_minus;
            K = P_minus*H'*inv(H*P_minus*H'+R_KF);
        end
        
        Sk = H*P_minus*H' + R_KF;
        x_hat_plus(:,k+1) = x_hat_minus(:,k+1) + K*innov;
        %Joseph formulation
        P_plus =(eye(4)-K*H)*P_minus*(eye(4)-K*H)' + K*R_KF*K';
        
        %compute NEES and NIS statistics
        NEESsshist(k) = (x_noisy(:,k)-x_hat_plus(:,k))'*inv(P_plus)*(x_noisy(:,k)-x_hat_plus(:,k));
        NISsshist(k) = innov'*inv(Sk)*innov / (length(innov)/3);
        
    end
    NEESamps(s,:) = NEESsshist;
    NISamps(s,:) = NISsshist;
    fprintf('s=%.0f \n',s)
end

%plot NEES statistics
epsNEESbar = mean(NEESamps,1);
alphaNEES = 0.05; %significance level
Nnx = Nsim*4; %N*n
r1x = chi2inv(alphaNEES/2,Nnx)./Nsim;
r2x = chi2inv(1-alphaNEES/2,Nnx)./Nsim;
figure; hold on; grid on; grid minor;
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2)
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES Statistics, avg \epsilon_x')
xlabel('time step k')
title(sprintf('NEES Estimation Results, N=%.0f',Nsim))
legend('NEES @ time k','r_1 bound','r_2 bound')
ylim([0 10])

epsNISbar = mean(NISamps,1);
alphaNIS = 0.05; %significance level
Nny = Nsim*3; %N*p
r1y = chi2inv(alphaNIS/2,Nny)./Nsim;
r2y = chi2inv(1-alphaNIS/2,Nny)./Nsim;
figure; hold on; grid on; grid minor;
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2)
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS Statistics, avg \epsilon_y')
xlabel('time step k')
title(sprintf('NIS Estimation Results, N=%.0f',Nsim))
legend('NIS @ time k','r_1 bound','r_2 bound')
ylim([0 10])

figure; hold on;
sgtitle('LKF Predicted States')
subplot(4,1,1); hold on; grid on; grid minor;
plot(tvec,x_hat_plus(1,:),'b-','LineWidth',2)
ylabel('X [km]')
subplot(4,1,2); hold on; grid on; grid minor;
plot(tvec,x_hat_plus(2,:),'b-','LineWidth',2)
ylabel('Xdot [km/s]')
subplot(4,1,3); hold on; grid on; grid minor;
plot(tvec,x_hat_plus(3,:),'b-','LineWidth',2)
ylabel('Y [km]')
subplot(4,1,4); hold on; grid on; grid minor;
plot(tvec,x_hat_plus(4,:),'b-','LineWidth',2)
ylabel('Ydot [km/s]')

%create ode45 simulation to compare against
s0 = x0 + dx0;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, x_perturbed] = ode45(@(t,s)orbit_prop_func(t,s),tvec,s0,opts);
x_perturbed=x_perturbed';

figure; hold on;
sgtitle('EKF State Estimation Errors')
subplot(4,1,1); hold on; grid on; grid minor;
plot(tvec,x_hat_plus(1,:)-x_perturbed(1,:),'b-','LineWidth',2)
ylabel('X [km]')
subplot(4,1,2); hold on; grid on; grid minor;
plot(tvec,x_hat_plus(2,:)-x_perturbed(2,:),'b-','LineWidth',2)
ylabel('Xdot [km/s]')
subplot(4,1,3); hold on; grid on; grid minor;
plot(tvec,x_hat_plus(3,:)-x_perturbed(3,:),'b-','LineWidth',2)
ylabel('Y [km]')
subplot(4,1,4); hold on; grid on; grid minor;
plot(tvec,x_hat_plus(4,:)-x_perturbed(4,:),'b-','LineWidth',2)
ylabel('Ydot [km/s]')










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

