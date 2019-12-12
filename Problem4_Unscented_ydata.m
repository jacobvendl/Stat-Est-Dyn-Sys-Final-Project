%Jake Vendl and Jack Toland

close all; clear all; clc

%Unscented KF

%load in provided data
load('orbitdeterm_finalproj_KFdata.mat')

%set problem inputs
mu = 398600;             % km^3/s^2
r0 = 6678;               % km
global rE wE
rE = 6378;               % km
wE = 2*pi/86400;         % rad/s
dt = 10;                 % s
P = 2*pi*sqrt(r0^3/mu);  % s
Q = eye(4)*1e-9; Q(1,1)=0; Q(3,3)=0;
R = eye(3)*1e-3; R(2,2)=.1;

%set x0
x0 = [6678, 0, 0, r0*sqrt(mu/r0^3)]';
dx0 = [0, 0.01, 0, 0.01]';

s0=x0+dx0;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, x_star] = ode45(@(t,s)orbit_prop_func(t,s),tvec,s0,opts);
x_star=x_star';

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, x_nom] = ode45(@(t,s)orbit_prop_func(t,s),tvec,x0,opts);
x_nom=x_nom';

x(:,1) = x0+dx0;
P_plus = eye(4)*1e6;

%set UKF inputs
L=4;  %states
m=3;  %measurements
alpha=1e-4;   %TUNE
ki=0;         %TUNE
beta=2;       %TUNE
lambda=alpha^2*(L+ki)-L;                    
c=L+lambda;                                
Wm=[lambda/c 0.5/c+zeros(1,2*L)];           
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);               
c=sqrt(c);
for k=1:length(tvec)-1
    %call function to create 9 sigma points
    X_hat=sigmas(x(:,k),P_plus,c);
    
    %call ukx to perform transform on state
    [x_hat_minus,X1,P_minus,devX]=utx(X_hat,Wm,Wc,L,Q,tvec,k);
    
    %call uky to perform transform on measurement
    [y_hat_minus,~,Pyy,devY]=uty(X1,Wm,Wc,m,R,tvec,k);
    
    if isempty(devY)==1
        x(:,k+1) = x_hat_minus;
        P_plus = P_minus;
    else 
        %calculate cross-covariance matrix
        Pxy=devX*diag(Wc)*devY';
        
        %5calculate Kalman gain
        K=Pxy*inv(Pyy);
        
        %find provided data measurement at time
        z = ydata{k+1}(1:3,1);
        s = size(ydata{k+1});
        if s(2) == 2 && length(y_hat_minus)>3 %two measurements recorded
            z=vertcat(z,ydata{k+1}(1:3,2));
        end
        if s(2) == 1 && length(y_hat_minus)>3
            y_hat_minus=y_hat_minus(1:3); %chop off measurements
            K = K(1:4,1:3);
            Pyy=Pyy(1:3,1:3);
        end
        
        %calculate state update
        x(:,k+1)=x_hat_minus+K*(z-y_hat_minus);
        %P_plus=P_minus-Pxy*inv(Pyy)*Pxy';          %covariance update
        P_plus = P_minus - K*Pyy*K';
        P_plus = (P_plus+P_plus')/2;
    end
    
    %save off 2 sigma values
    twoSigX(k+1) = 2*sqrt(P_plus(1,1));
    twoSigXdot(k+1) = 2*sqrt(P_plus(2,2));
    twoSigY(k+1) = 2*sqrt(P_plus(3,3));
    twoSigYdot(k+1) = 2*sqrt(P_plus(4,4));
    fprintf('k=%.0f\n',k)
end

%make plots
figure; hold on; 
sgtitle('UKF Predicted States')
subplot(4,1,1); hold on; grid on; grid minor;
plot(tvec,x(1,:),'b-','LineWidth',2)
ylabel('X [km]')
subplot(4,1,2); hold on; grid on; grid minor;
plot(tvec,x(2,:),'b-','LineWidth',2)
ylabel('Xdot [km/s]')
ylim([-10 10])
subplot(4,1,3); hold on; grid on; grid minor;
plot(tvec,x(3,:),'b-','LineWidth',2)
ylabel('Y [km]')
subplot(4,1,4); hold on; grid on; grid minor;
plot(tvec,x(4,:),'b-','LineWidth',2)
ylabel('Ydot [km/s]'); xlabel('Time [s]')
ylim([-10 10])

%create ode45 simulation to compare against
s0 = x0 + dx0;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, x_perturbed] = ode45(@(t,s)orbit_prop_func(t,s),tvec,s0,opts);
x_perturbed=x_perturbed';

figure; hold on;
sgtitle('UKF State Estimation Errors')
subplot(4,1,1); hold on; grid on; grid minor;
plot(tvec,x(1,:)-x_perturbed(1,:),'b-','LineWidth',2)
ylabel('X [km]')
subplot(4,1,2); hold on; grid on; grid minor;
plot(tvec,x(2,:)-x_perturbed(2,:),'b-','LineWidth',2)
ylabel('Xdot [km/s]')
ylim([-1.5 1.5])
subplot(4,1,3); hold on; grid on; grid minor;
plot(tvec,x(3,:)-x_perturbed(3,:),'b-','LineWidth',2)
ylabel('Y [km]')
subplot(4,1,4); hold on; grid on; grid minor;
plot(tvec,x(4,:)-x_perturbed(4,:),'b-','LineWidth',2)
ylabel('Ydot [km/s]'); xlabel('Time [s]')
ylim([-1.5 1.5])


function [x1,X1,P1,X2]=utx(X,Wm,Wc,n,Q,tvec,t)
%ut for state
L=size(X,2);
x1=zeros(n,1);      %set transformed average
X1=zeros(n,L);      %set transformed samples
for k=1:L
    %deterministic function for finding x
    [~,temp] = ode45(@(t,s)orbit_prop_func(t,s),[tvec(t) tvec(t+1)],X(:,k),odeset('RelTol',1e-12,'AbsTol',1e-12));
    X1(:,k)=temp(end,:);
    x1=x1+Wm(k)*X1(:,k);
end
X2=X1-x1(:,ones(1,L));  %set transformed deviations
P1=X2*diag(Wc)*X2'+Q;   %set transformed covariance
end

function [z1,Z1,P2,Z2]=uty(X1,Wm,Wc,m,R,tvec,t)
global rE wE
%ut for measurements
L=size(X1,2);

z1=[];      %set transformed average
Z1=[];      %set transformed samples
for k=1:L
    %deterministic function for finding y
    [~,temp] = ode45(@(t,s)orbit_prop_func(t,s),[tvec(t) tvec(t+1)],X1(:,k),odeset('RelTol',1e-12,'AbsTol',1e-12));
    %use current state estimate to predict measurement perfectly at time t
    X=temp(end,1); XD=temp(end,2); Y=temp(end,3); YD=temp(end,4);
    meas=[];
    for i=1:12
        theta = (i-1)*pi/6;
        currentTime = tvec(t);
        
        %find station position and velocity
        Xs = rE*cos(wE*currentTime + theta);
        Ys = rE*sin(wE*currentTime + theta);
        XDs = -rE*wE*sin(wE*currentTime + theta);
        YDs = rE*wE*cos(wE*currentTime + theta);
        
        %process to see if a given station can see the s/c
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
            
            rho = sqrt((X-Xs)^2 + (Y-Ys)^2);
            rhoDot = ((X-Xs)*(XD-XDs) + (Y-Ys)*(YD-YDs)) / rho;
            meas = vertcat(meas,[rho;rhoDot;phi]);
        end
    end
    if isempty(meas)==1
        Z1=[];
        z1=[];
        break
    end
    Z1(:,k)=meas;
    if k==1
        z1 = zeros(length(meas),1);
    end
    z1=z1+Wm(k)*Z1(:,k);
end
%perform checks to make sure a measurement came through and handle the case
%where it doesn't
if isempty(Z1)==1
    Z2=[];
    P2=[];          
else
    Z2=Z1-z1(:,ones(1,L));          %set transformed deviations
    if length(meas)>3
        R = blkdiag(R,R);
    end
    P2=Z2*diag(Wc)*Z2'+R;           %set transformed covariance
end
end

function X=sigmas(x,P,c)
%calculate sigma points around the reference point
%NOTE: c=sqrt(L+lambda), that step has already been handled
A = c*chol(P)'; %perform chol decomp
Y = x(:,ones(1,numel(x)));  %lay framework for sigma distribution based on given reference point
X = [x Y+A Y-A];    %calculate 9 sigma points based on chol decomposition
end

function [ ds ] = orbit_prop_func(t,s)
%propagation function for motion
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
