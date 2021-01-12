%% Load data
load('..\Assignment_Data_SC42145.mat')
%% Define the uncertainty block
%ASK
blk=[ 1 1; 1 1; 1 1; 1 1];
%% Part 1
s=tf('s');

% Define Uncertainty Weights
Wi1=((s/16*pi)+0.3)/((s/64*pi)+1);
Wi2=((s/16*pi)+0.3)/((s/64*pi)+1);
Wo1=(0.05*s+0.2)/(0.01*s+1);
Wo2=(0.05*s+0.2)/(0.01*s+1);

Wo=[Wo1 0; 0 Wo2];
Wi=[Wi1 0; 0 Wi2];
%Define Performance Weigths
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

%% Part 2

%% Part 3


% Bode Plots
figure(1)
bodemag(Wi1);
figure(2)
bodemag(Wi2);

% Maybe caused by unmodelled hihg-frequency dynamics (pg 269)

G = FWT(1:2,1:2);
unc = ultidyn('unc',[1 1]); 

Di=[ultidyn('Di1',[1,1]) 0; 0 ultidyn('Di2',[1,1])];
Do=[ultidyn('Do1',[1,1]) 0; 0 ultidyn('Do2',[1,1])];

usys = (eye(2) + Wo*Do)*G*(eye(2) + Wi*Di);

%Plot of the singular values
figure(3)
sigma(usys)

%Plot of comparison Bode-Nominal Plant
figure(4)
bodemag(usys,G) %Change linewidth of the nominal

%% Generalized Plant

% Construction of the generalized plant
systemnames ='G Wp Wu Wi Wo'; % Define systems
inputvar ='[u_d(2);u_od(2);w(2);u(2)]'; % Input generalized plant
input_to_G= '[u_d+u]';
input_to_Wu= '[u]';
input_to_Wp= '[-w+G+u_od]';
input_to_Wi= '[u]';
input_to_Wo= '[G]';
outputvar= '[Wi; Wo;Wp; Wu; w-G]'; % Output generalized plant
sysoutname='P'; cleanupsysic= 'yes'; 
sysic;
P=minreal(P);

% Pm=[zeros(2) zeros(2) zeros(2) Wi; Wo*G zeros(2) zeros(2)  Wo*G;...
%     Wp*G Wp -Wp Wp*G; zeros(2) zeros(2) zeros(2) Wu; -G -eye(2) eye(2) -G];

%% Part 5
%Load Controller designed in PART2
load('ControllerPart2.mat')

% Remember to do comments on mu!!!!!!!!
N=lft(P,K);
max(real(eig(N))) %Since the larger real part is negative, N is NS.

% Open loop transfer matrix
L = minreal(G*K);

% Generalized Nyquist plot
IpL = minreal(eye(2)+L);
det_IpL = IpL(1,1)*IpL(2,2)-IpL(2,1)*IpL(1,2);
Lpoles = pole(L);
figure()
nyquist(minreal(det_IpL));

%%
%Define frequency range
omega=logspace(-5,3,50*61);

% Define Nf
Nf=frd(N,omega);

% Check for NP

N22=N(7:8,5:6);
Dp=ultidyn('Dp',[2,4]);

% N_NP=lft(Dp,N22);
% sigma(N_NP)
norm(N22,'inf') %The maximum value of the inf norm is 0.6843, which...
...coincides with the largest singular value



blk=[ 2 4]; % Full complex uncertainty block
[mubnds,muinfo]=mussv(Nf(5:8,5:6),blk,'c');
muNP=mubnds(:,1);
[muNPinf, muNPw]=norm(muNP,inf);
    
% Check for RS

%Create a combined uncertainty matrix



blk=[ 1 1; 1 1;  1 1; 1 1]; % structured uncertainty
[mubnds,muinfo]=mussv(Nf(1:4,1:4),blk,'c');
muRS=mubnds(:,1);
[muRSinf, muRSw]=norm(muRS,inf);

 %YEZ RS :D
%%
% Check for RP



blk=[ 1 1; 1 1; 1 1; 1 1; 2 4];
% blk=[ 1 1; 1 1; 1 1; 1 1; 2 2];
[mubnds,muinfo]=mussv(Nf,blk,'c');
muRP=mubnds(:,1);
[muRPinf, muRPw]=norm(muRP,inf);
 %Nope RP :(
 
 %% Part 3.2
 
 % Part 1 : haha no
 
 %% Part 2: D-K Iterations
 
 % D-K iterations auto-tuning
% Delta=[ultidyn('D 1',[1,1]) 0; 0 ultidyn('D 2',[1,1])];
D=[Di zeros(2);zeros(2) Do];
Punc=lft(D,P);
% opt=dkitopt('FrequencyVector', omega,'DisplayWhileAutoIter','on')
opt=dkitopt('DisplayWhileAutoIter','on');
[Ki1,clp,dkinfo]=musyn(Punc,2,2,K);

%% Manuel

systemnames ='G Wp Wu Wi Wo'; % Define systems
inputvar ='[u_d(2);u_od(2);w(2);u(2)]'; % Input generalized plant
input_to_G= '[u_d+u]';
input_to_Wu= '[u]';
input_to_Wp= '[-w+G+u_od]';
input_to_Wi= '[u]';
input_to_Wo= '[G]';
outputvar= '[Wi; Wo;Wp; w-G]'; % Output generalized plant
sysoutname='P'; cleanupsysic= 'yes'; 
sysic;
P=minreal(P);

 blk=[ 1 1; 1 1;  1 1; 1 1; 2 2];
[K2,CL,GAM,INFO] = hinfsyn(P,2,2);

for i=1:1:10
Nf=frd(lft(P,K2),omega);
[mubnds,muinfo]=mussv(Nf,blk,'c');
muRP=mubnds(:,1); [muRPinf, muRPw]=norm(muRP,inf);
[VDelta,VSigma,VLmi] = mussvextract(muinfo);
D=VSigma.DLeft;
dd1 = fitmagfrd((D(1,1)/D(3,3)),6);
dd2 = fitmagfrd((D(2,2)/D(3,3)),6);
dd3 = fitmagfrd((D(4,4)/D(6,6)),6);
dd4 = fitmagfrd((D(5,5)/D(6,6)),6);
Dscale=minreal(append(dd1, dd2,dd3, dd4,tf(eye(4))));
[K2,CL,GAM2,INFO] = hinfsyn(minreal(Dscale*P*inv(Dscale)),2,2);
end

%% Time-domain Simulation
% Extends the original system with a direct feedthrough for the V
Kex=[K2 zeros(2,1);0 0 1];
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
