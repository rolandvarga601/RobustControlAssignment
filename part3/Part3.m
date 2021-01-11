%% Load data
myfile = 'Assignment_Data_SC42145.mat';
[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,myfile))
%% Define the uncertainty block

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

sysNom = FWT(1:2,1:2);
unc = ultidyn('unc',[1 1]); 

Delta_i=[ultidyn('Di1',[1,1]) 0; 0 ultidyn('Di2',[1,1])];
Delta_o=[ultidyn('Do1',[1,1]) 0; 0 ultidyn('Do2',[1,1])];

usys = (1 + Wo*Delta_o)*sysNom*(1 + Wi*Delta_i);



% % Set properties of usys
% usys.InputName = 'u';
% usys.OutputName = 'fs';

