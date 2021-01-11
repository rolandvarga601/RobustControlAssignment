%% Init
clear all
load Assignment_Data_SC42145.mat

%% Defining the MIMO system
g_mimo = tf(FWT(1:2,1:2));

%% Calculating the RGA
g_mimo_freq1=evalfr(g_mimo,0*1i); %Evaluate g_mimo at w=0
RGA_freq1=(g_mimo_freq1).*(inv(g_mimo_freq1))'; %Calculate RGA at w=0
g_mimo_freq2=evalfr(g_mimo,0.4*2*pi*1i); %Evaluate g_mimo at w=0.4*2*pi
RGA_freq2=(g_mimo_freq2).*(inv(g_mimo_freq2))'; %Calculate RGA at w=0.4*2*pi

%% Calculating the MIMO poles and zeros
Z1=tzero(minreal(balreal(g_mimo)));
P1=pole(minreal(balreal(g_mimo)));

%% Designing the weights
s=tf('s');

% Output error weight
wb=0.4*2*pi; 
A=10^-4;
M=1.8;
Wp11=((s/M)+wb)/(s+wb*A);
Wp=[Wp11 0; 0 0.2];
%bode(Wp11)

% Control input weight
Wu=[0.01 0;...
    0 (5*(10^-3)*s^2+7*(10^-4)*s+5*(10^-5))/(s^2+14*10^-4*s+10^-6)];

%% General plant construction and Hinf synthesis
% Changing notation
G=g_mimo;

% Construction of the generalized plant
systemnames ='G Wp Wu '; % Define systems
inputvar ='[w(2);u(2)]'; % Input generalized plant
input_to_G= '[u]';
input_to_Wu= '[u]';
input_to_Wp= '[w-G]';
outputvar= '[Wp; Wu; w-G]'; % Output generalized plant
sysoutname='P';
sysic;

% Hinf controller synthesis
ncont = 2; 
nmeas = 2;
[K,CL,GAM,INFO] = hinfsyn(P,nmeas,ncont);
K = minreal(K);

%% System analysis
% Open loop transfer matrix
L = minreal(G*K);

% Generalized Nyquist plot
IpL = minreal(eye(2)+L);
det_IpL = IpL(1,1)*IpL(2,2)-IpL(2,1)*IpL(1,2);
Lpoles = pole(L);
figure()
nyquist(minreal(det_IpL));

% Internal stability check
Q = minreal(K*inv(eye(2)+G*K));
det_Q = minreal(Q(1,1)*Q(2,2)-Q(2,1)*Q(1,2));
figure()
nyquist(det_Q);

%% Simulation - reference tracking
% Extends the original system with a direct feedthrough for the V
Kex=[K zeros(2,1);0 0 1];
sys_fb = feedback(FWT*Kex,eye(2),[1 2],[1 2]);

% System response to step in omega reference
figure()
step(sys_fb(1,1))

% Transfer matrix from reference input to control signal output
sys_fb_control = feedback(Kex,FWT,[1 2],[1 2 3]);

% Control inputs in case of step in omega reference
figure()
step(sys_fb_control(1,1))
figure()
step(sys_fb_control(2,1))

%% Simulation - disturbance rejection
% System response to step in the disturbance
figure()
step(sys_fb(1,3))

% Control input response to step in the disturbance
figure()
step(sys_fb_control(1,3))
figure()
step(sys_fb_control(2,3))

%% Part 2.1
% Wp weight design
wb=10*2*pi;
A=10^-4;
M=1.8;
Wp11=1/M*(s+wb*M)/(s+wb*A);
Wpd = Wp11;

% Wu weight design
A=10^-4;
M=1.8;
% wbc1=0.005*2*pi;      % Initial design
% wbc2=0.01*2*pi;       % Initial design
wbc1=0.002*2*pi;        % Final design
wbc2=0.05*2*pi;         % Final design
Wud=[(1/(M*10)*(s/(wbc1*M/10)+1)/(s/(wbc1/A)+1)) 0;...
            0 (((s/(M*20))+wbc2)/(s+wbc2*A))];

% Redefining G for simplicity
G = G(1,:);

% Creating the SISO disturbance transfer function
Gd = minreal(tf(FWT(1,3)));

% Contruction of the generalized plant
systemnames ='G Gd Wpd Wud '; % Define systems
inputvar ='[w(1);u(2)]'; % Input generalized plant
input_to_G= '[u]';
input_to_Gd= '[w]';
input_to_Wud= '[u]';
input_to_Wpd= '[-Gd-G]';
outputvar= '[Wpd; Wud; -Gd-G]'; % Output generalized plant
sysoutname='P';
sysic;

% Controller synthesis
[Kd,CLd,GAMd,INFOd] = hinfsyn(P,1,2);
Kd = minreal(Kd);

% Construction of the closed loop system
CL_dist = minreal(feedback(1,G*Kd)*Gd);     %==S*Gd
CL_control = minreal(feedback(Kd,G)*Gd);    %==Kd*S*Gd

%% Comparing the weights to the closed loop transfer functions
% Performance comparison
figure()
bode(1/Wpd,CL_dist)

% Control input comparison
figure()
bode(1/Wud(1,1),CL_control(1,1),1/Wud(2,2),CL_control(2,1))

% feedback(Kd,G) <=> Kd*feedback(eye(2),G*Kd) <=> Kd*S

%% Simulation - real-life data
% Interpolation of the given disturbance signal
xi = linspace(0,max(Wind_Data.time),size(Wind_Data.time,1));
yi = interp1(Wind_Data.time,Wind_Data.data,xi,'linear');

% System response to the given disturbance
figure()
[omega_real,t]=lsim(CL_dist, yi, xi);
plot(t,omega_real)

% Control input response to the given disturbance
[input_real,t]=lsim(CL_control, yi, xi);
figure()
plot(t,input_real(:,1))
figure()
plot(t,input_real(:,2))

%% Simulation - low frequency disturbance
% Creating "filtered" signal
t = linspace(0,max(Wind_Data.time),size(Wind_Data.time,1)*10);
Wind_Data_filtered = 2*sin(t*2*pi/1000);    % A=2[-] T=1000[s]

% Comparison to see if the filtered signal fit the original
%figure()
%plot(xi,Wind_Data.data,t,Wind_Data_filtered)

% System response to the disturbance
figure()
[omega_noiseless,t]=lsim(CL_dist, Wind_Data_filtered, t);
plot(xi,omega_real,t,omega_noiseless)

CL_control = minreal(Gd*feedback(Kd,G));    %==Gd*Kd*S

% System response to the disturbance
[input_noiseless,t]=lsim(CL_control, Wind_Data_filtered, t);
figure()
plot(xi,input_real(:,1),t,input_noiseless(:,1))
figure()
plot(xi,input_real(:,2),t,input_noiseless(:,2))
