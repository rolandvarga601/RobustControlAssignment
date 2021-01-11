%% Defining the MIMO system
clear all
load Assignment_Data_SC42145.mat

%% Defining the mimo system
g_mimo = tf(FWT(1:2,1:2));

%% Calculating the RGA
g_mimo_freq1=evalfr(g_mimo,0*1i); %Evaluate g_mimo at w=0
RGA_freq1=(g_mimo_freq1).*(inv(g_mimo_freq1))'; %Calculate RGA at w=0
g_mimo_freq2=evalfr(g_mimo,0.4*2*pi*1i); %Evaluate g_mimo at w=0.4*2*pi
RGA_freq2=(g_mimo_freq2).*(inv(g_mimo_freq2))'; %Calculate RGA at w=0.4*2*pi

%% Calculating the MIMO poles and zeros
Z1=tzero(minreal(balreal(g_mimo)));
P1=pole(minreal(balreal(g_mimo)));

%% Designing Wp11
s=tf('s');
wb=0.4*2*pi; 
A=10^-4;
M=1.8;
Wp11=((s/M)+wb)/(s+wb*A);
Wp=[Wp11 0; 0 0.2];
%bode(Wp11)
%% Part 7
Wu=[0.01 0; 0 (5*(10^-3)*s^2+7*(10^-4)*s+5*(10^-5))/(s^2+14*10^-4*s+10^-6)];
P_daniel=[Wp Wp*g_mimo; zeros(2) -Wu; -eye(2) -g_mimo];
P_compare=[Wp*g_mimo; -Wu; -g_mimo];
% P_ss=tf2ss(P_compare);
G=g_mimo;
P_aug=augw(G,Wp,Wu,[]);
P_aug=tf(augw(G,Wp,Wu,[]));
P_aug(:,1:2)=0*s;
ncont = 2; 
nmeas = 2; 
[K1,CL1,gamma,info] = hinfsyn(P_aug,nmeas,ncont);



systemnames ='G Wp Wu '; % Define systems
inputvar ='[w(2);u(2)]'; % Input generalized plant
input_to_G= '[u]';
input_to_Wu= '[u]';
input_to_Wp= '[w-G]';
outputvar= '[Wp; Wu; w-G]'; % Output generalized plant
sysoutname='P';
sysic;
[K,CL,GAM,INFO] = hinfsyn(P,2,2);
K = minreal(K);
% Wu2=tf([0.01 0; 0 0.00001]);
% systemnames ='G Wp Wu2 '; % Define systems
% inputvar ='[w(2);u(2)]'; % Input generalized plant
% input_to_G= '[u]';
% input_to_Wu2= '[u]';
% input_to_Wp= '[w-G]';
% outputvar= '[Wp; Wu2; w-G]'; % Output generalized plant
% sysoutname='P';
% sysic;
% [K,CL,GAM,INFO] = hinfsyn(P,2,2);

%[K2,CL2,GAM2,INFO2]=mixsyn(G,Wp,Wu,[]);

%% System analysis
L = minreal(G*K);
IpL = minreal(eye(2)+L);
det_IpL = IpL(1,1)*IpL(2,2)-IpL(2,1)*IpL(1,2);
Lpoles = pole(L);
figure()
nyquist(minreal(det_IpL));

% Internal stability check L3S6
Q = minreal(K*inv(eye(2)+G*K));
det_Q = minreal(Q(1,1)*Q(2,2)-Q(2,1)*Q(1,2));
figure()
nyquist(det_Q);

%% Simulation - reference tracking
% Extends the original system with a direct feedthrough for the V
Kex=[K zeros(2,1);0 0 1];
sys_fb = feedback(FWT*Kex,eye(2),[1 2],[1 2]);
figure()
step(sys_fb(1,1))

sys_fb_control = feedback(Kex,FWT,[1 2],[1 2 3]);
figure()
step(sys_fb_control(1,1))
figure()
step(sys_fb_control(2,1))

%% Simulation - disturbance rejection
figure()
step(sys_fb(1,3))

figure()
step(sys_fb_control(1,3))
figure()
step(sys_fb_control(2,3))

%% Part 2.1

% Weight design
%wb=1*2*pi;
wb=10*2*pi;
A=10^-4;
M=1.8;
Wp11=1/M*(s+wb*M)/(s+wb*A);
Wpd = Wp11;
%Wud = [1/Wp11 0; 0 Wp11];
A=10^-4;
M=1.8;
% wbc1=0.005*2*pi;
% wbc2=0.01*2*pi;
% wbc1=0.002*2*pi;
% wbc2=0.01*2*pi;
% wbc2=0.5*2*pi;
% wbc1=0.001*2*pi;
wbc1=0.002*2*pi;
wbc2=0.05*2*pi;
%Wud=[(1/M*(s/(wbc1*M)+1)/(s/(wbc1/A)+1)) 0; 0 (((s/M)+wbc2)/(s+wbc2*A))];
%Wud=[(1/(M*10)*(s/(wbc1*M/10)+1)/(s/(wbc1/A)+1)) 0; 0 (((s/(M*10))+wbc2)/(s+wbc2*A))];
%Wud=[(1/(M*10)*(s/(wbc1*M/10)+1)/(s/(wbc1/A)+1)) 0; 0 (((s/(M*20))+wbc2)/(s+wbc2*A))];
Wud=[(1/(M*10)*(s/(wbc1*M/10)+1)/(s/(wbc1/A)+1)) 0; 0 (((s/(M*20))+wbc2)/(s+wbc2*A))];
G = G(1,:);
Gd = minreal(tf(FWT(1,3)));
systemnames ='G Gd Wpd Wud '; % Define systems
inputvar ='[w(1);u(2)]'; % Input generalized plant
input_to_G= '[u]';
input_to_Gd= '[w]';
input_to_Wud= '[u]';
input_to_Wpd= '[-Gd-G]';
outputvar= '[Wpd; Wud; -Gd-G]'; % Output generalized plant
sysoutname='P';
sysic;
[Kd,CLd,GAMd,INFOd] = hinfsyn(P,1,2);
Kd = minreal(Kd);

%% Simulation
CL_dist = minreal(feedback(1,G*Kd)*Gd); %==S*Gd
xi = linspace(0,max(Wind_Data.time),size(Wind_Data.time,1));
yi = interp1(Wind_Data.time,Wind_Data.data,xi,'linear');
figure()
% %step(CL_dist)
[omega_real,t]=lsim(CL_dist, yi, xi);
plot(t,omega_real)

CL_control = minreal(feedback(Kd,G)*Gd);    %==Gd*Kd*S
[input_real,t]=lsim(CL_control, yi, xi);
figure()
plot(t,input_real(:,1))
figure()
plot(t,input_real(:,2))

%% Bode plots
figure()
bode(1/Wpd,CL_dist)
figure()
bode(1/Wud(1,1),CL_control(1,1),1/Wud(2,2),CL_control(2,1))

% feedback(Kd,G) <=> Kd*feedback(eye(2),G*Kd) <=> Kd*S
%%
%t = Wind_Data.time;
t = linspace(0,max(Wind_Data.time),size(Wind_Data.time,1)*10);
Wind_Data_filtered = 2*sin(t*2*pi/1000);    % A=2[-] T=1000[s]
%Wind_Data_filtered = 2*sin(xi*2*pi/1000);%*pi);
%plot(t,Wind_Data.data,t,Wind_Data_filtered)

figure()
% %step(CL_dist)
[omega_noiseless,t]=lsim(CL_dist, Wind_Data_filtered, t);
plot(xi,omega_real,t,omega_noiseless)

CL_control = minreal(Gd*feedback(Kd,G));    %==Gd*Kd*S
[input_noiseless,t]=lsim(CL_control, Wind_Data_filtered, t);
figure()
plot(xi,input_real(:,1),t,input_noiseless(:,1))
figure()
plot(xi,input_real(:,2),t,input_noiseless(:,2))




