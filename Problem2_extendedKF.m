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

for k=1:length(tvec)
    
end