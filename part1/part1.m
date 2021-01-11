%% Defining the SISO system
clear all
load Assignment_Data_SC42145.mat

g_siso = tf(FWT(1,1));

% Define sys as -1 times the original system (because of phase problems)
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

%% No need for this but just in case
%[mag,phase,wout] = bodeplot(sys,opts);
% mag = squeeze(mag);
% phase = squeeze(phase);
% figure()
% loglog(wout,mag)
% figure()
% semilogx(wout,phase-ones(size(phase))*360)

%% PI controller design
opts = bodeoptions('cstprefs'); % Load System Id Toolbox default
%opts.FreqUnits = 'Hz';
opts.PhaseMatching = 'on'; % Bode plot will start at 0 deg instead of 360

s = tf('s');
Kpi = 1;
%K =tf(0.90903*(s+0.30995)/s); % PI tuner controller
K =tf(Kpi*(s+0.203)/s);
figure('Name','Open loop bode','NumberTitle','off')
bodeplot(K*sys,opts)

%% Choosing the PI controller gain to get 60 deg PM (120 deg phase)
Kpi = 1/0.192;
figure('Name','60 deg PM','NumberTitle','off')
margin(Kpi*K*sys);
figure('Name','Step response -- 60 deg PM','NumberTitle','off')
step(feedback(Kpi*K*sys,1)) % Too noisy

%% Choosing the PI controller gain to get 95 deg PM (85 deg phase)
%Kpi = 10^(2.3/20); % Works good
Kpi = 1/0.309;
K =tf(Kpi*(s+0.203)/s);
figure('Name','95 deg PM','NumberTitle','off')
margin(K*sys);
figure('Name','Step response -- 95 deg PM','NumberTitle','off')
step(feedback(Kpi*K*sys,1)) % Still too much overshoot
% 0.03 Hz component in not supressed enough -> Increase gain margin

%% Increasing GM till overshoot is acceptable
Kpi = 10^(2.3/20); % Works good
K =tf(Kpi*(s+0.203)/s);
figure('Name','17.1 dB GM','NumberTitle','off')
margin(K*sys);
figure('Name','Step response -- 17.1 dB GM','NumberTitle','off')
step(feedback(K*sys,1))

%% Invert the controller sign
K = -K;
figure()
step(feedback(K*g_siso,1))

%% Complementary sensitivity function

%% Disturbance transfer function
% A_siso_dist = A;
% B_siso_dist = B(:,3);
% C_siso_dist = C(1,:);
% D_siso_dist = D(1,3);

feedin = [1];
feedout = [1];

% Converting to transfer function
% [b,a] = ss2tf(A_siso_dist,B_siso_dist,C_siso_dist,D_siso_dist);
Gd = tf(FWT(1,3));

% Closed loop system with the controller
Gd_cl = Gd*feedback(1,K*g_siso);

%FWT_cl = feedback(FWT,K,feedin,feedout,-1);
step(Gd_cl)
%step(FWT_cl(1,3))
