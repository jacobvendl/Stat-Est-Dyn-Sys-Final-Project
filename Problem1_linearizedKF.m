%Jake Vendl and Jack Toland
%ASEN 5044 Final Project

clear all; close all; clc

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
P = 2*pi*sqrt(r0^3/mu);  %s

%CODE OUTLINE: jack please check and maybe make your own? this is critical
%monte carlo to pick a dx
%simulate ground truth state using given Q
%simulate measurements using R
%run KF on the dx, saving out state estimation errors
%NEES and NIS tests

%STEP ONE  - generate input to truth model
x0 = [6678, 0, 0, r0*sqrt(mu/r0^3)];
dx0 = [0.1, 0.001, 0.1, 0.001]; %come back and get this with MC

%STEP TWO - simulate ground truth state using ode45
s0 = x0 + dx0;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[T, S] = ode45(@(t,s)orbit_prop_func(t,s),tvec,s0,opts);

%question stopping me: how do we use Q, the process noise?
    %my first idea is to tack on w(k) for each time step









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
