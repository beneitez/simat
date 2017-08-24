close all
clear all
clc

% Data analysis for the amplitude file
% 20160727 - Miguel Beneitez - beneitez@kth.se

writefile = 1;

% Fixed parameters for the simulation
Retau  = 180;
Recl   = 4200;

% Eddy Turnover Time for the simulation. Non dimensional
% lifetime of an eddy, based on the half channel and the 
% nondimensional u_tau
utau    = Retau/Recl;
h       = 1;
ett     = h/utau;

% Data reading from the amplitude file

matAmp = load('amp.stat');

t      = matAmp(:,1);
urms   = matAmp(:,2);
vrms   = matAmp(:,3);
wrms   = matAmp(:,4);
voxrms = matAmp(:,5);
voyrms = matAmp(:,6);
vozrms = matAmp(:,7);
dUuv   = matAmp(:,9);
e0     = matAmp(:,10);
hplus  = matAmp(:,11);
px     = matAmp(:,12);
Reb    = matAmp(:,13);
tau1   = matAmp(:,14);
tau2   = matAmp(:,15);

% Locating the position for which we start control
% Simulations carried out for start at t=500

[~,I] = min(abs(t-500));

if I~=1
   t      = t(I:end);
   urms   = urms(I:end);
   vrms   = vrms(I:end);
   wrms   = wrms(I:end);
   voxrms = voxrms(I:end);
   voyrms = voyrms(I:end);
   vozrms = vozrms(I:end);
   dUuv   = dUuv(I:end);
   e0     = e0(I:end);
   hplus  = hplus(I:end);
   px     = px(I:end);
   Reb    = Reb(I:end);
   tau1   = tau1(I:end);
   tau2   = tau2(I:end);
end

% Time shift to have developed TCF at t=0

t  = t-500;
[~,Ie] = min(abs(t/ett-31));

plot(t(1:Ie)/ett,-px(1:Ie))
axis([t(1)/ett 30 1.4e-3 2e-3])

if writefile==1
    tvsPx = fopen('tVsPx.txt','w','n');
    fprintf(tvsPx,'%2.12f %2.12f \n',[t(1:Ie)/ett -px(1:Ie)]');
end