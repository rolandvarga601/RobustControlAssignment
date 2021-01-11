%% Defining the SISO system
clear all
load Assignment_Data_SC42145.mat

% Extracting the transfer function between omega and beta
g_siso = tf(FWT(1,1));

% Define sys as -1 times the original system (because of convinience)
sys = zpk(-g_siso);
sys.DisplayFormat='frequency';

%% Analysis of the system
% Modify the default bodeplot options
opts = bodeoptions('cstprefs'); % Load System Id Toolbox default
opts.FreqUnits = 'Hz';
opts.PhaseMatching = 'on'; % Bode plot will start at 0 deg instead of 360

figure('Name','System Bode plot','NumberTitle','off')
bodeplot(sys,opts);

figure('Name','System pole-zero map','NumberTitle','off')
pzmap(sys);

%% I controller design
opts = bodeoptions('cstprefs'); % Load System Id Toolbox default
%opts.FreqUnits = 'Hz';
opts.PhaseMatching = 'on'; % Bode plot will start at 0 deg instead of 360

s = tf('s');
K =tf((0.203)/s);
figure('Name','Open loop bode','NumberTitle','off')
bodeplot(K*sys,opts)

%% Decreasing the GM till overshoot is acceptable and the system is fast
K =tf((0.26)/s); % Works good
figure('Name','16.5 dB GM','NumberTitle','off')
margin(K*sys);
figure('Name','Step response -- 16.5 dB GM','NumberTitle','off')
step(feedback(K*sys,1))

% Printing the time-domain informations
stepinfo(feedback(K*sys,1))

%% Inverting the controller sign
K = -K;
figure()
step(feedback(K*g_siso,1))

%% Exercise 4: Disturbance rejection test
% Extracting Gd from the original system
Gd = tf(FWT(1,3));

% Closed loop system with the controller
Gd_cl = Gd*feedback(1,K*g_siso);

figure('Name','Disturbance rejection','NumberTitle','off')
step(Gd_cl)

