% apply damping to irregular wave S2, Fwd, A30

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

%% load experimental data
dis = [35 22 9]*1e-3;
phi = [10 30 50];
direction = {'Fwd', 'Aft'};                 % wave direction [0, pi]
Tp = [1 1.25 1.5 1 1.25 1.5];               % peak wave period
Hs = [10.6 18.3 29.6 16.9 28.4 45.3]/1000;  % significant wave height

myCell=cell(18,1); % cell to hold plot labels
j=0;

k = 1; % forward direction
m = 2; %:length(Tp)
n=1; %:length(phi)
       
%inspect

%% load experimental data
dataPath = '/home/s2122199/Documents/Edinburgh/mocean/data/UoE Flume Tests 2019-02/IrrTS';
modelledDataPath = '/home/s2122199/Documents/Edinburgh/mocean/data/modelled_irregular_data/';
dataSetName = [direction{k} '_A' num2str(phi(n)) '_Irr00TS' num2str(m)];
forceTab = readtable([dataPath '/' dataSetName '.csv']);
            
expForces(:,1) = forceTab.Var1; % surge
expForces(:,2) = forceTab.Var3; % heave
expForces(:,3) = forceTab.Var5; % pitch
            
%% take fft of experimental timeseries
            
df = 64;
time = (0:length(expForces(:,3)) -1)/df; % experiment run time
timeSeries = expForces; % time series' of forces
tda = TimeDomainAnalysis;
tda.SetMotions(1, time, (timeSeries'));           % sets Motions as array Ntest x Ndof X Ntime. Input as DOF indices, time vector, data
[specs, freqs] = tda.GetMotions([], [], 'spectra', 'noMean', 'smooth', 0.1);
            
            
%% make power spectrum
            
% frequency from 0.5 Hz to 1.5 Hz vector to line up with fr
freqDat = freqs{3};
step = 12/time(end); % increase step size for comparison with numerical
Fstart = freqDat(round(time(end)*0.4));
Fend = freqDat(round(time(end)*1.8));
F0 = Fstart:step:Fend;
T0=1./flip(F0);
inds = ismembertol(freqDat,F0,1E-5); % get the indices of the experimental data to compare
            
            
expSpec =zeros(sum(inds),3);
dfD = freqDat(2) - freqDat(1);
for ff = 1:3
    FexDden = 0.5*abs(specs{ff}).^2./dfD; % force energy density spectrum
    expSpec(:,ff)= FexDden(inds);
end
FexExp(n+j,:,:) = expSpec;
            
myCell{n+j} = ['TS' num2str(m) ' ' num2str(phi(n)) char(176) ];
            


%% inspect the undamped modelled data
load([modelledDataPath dataSetName '_modelled']);
            
% get the undamped forces
FexU = squeeze(CompU.Fex); 
force_modelled=abs(FexU);
% make synthetic wave spectrum
spectra = Bretschneider(Hs(m), Tp(m), T0); % wave spectrum
A = spectra.Amplitudes; % amplitude spectum
            
Fs = flip(FexU.*A'); % force amplitude spectum
df = F0(2) - F0(1);
FD = 0.5*abs(Fs).^2./df; % force energy density spectrum
FexMod(n+j,:,:) = FD;
% surge force: DOF = 1
% heave force: DOF = 3
% pitch moment: DOF = 5


%% apply the damping plate method
%apply damping 

coef = ViscDampCoef;
 
coef.Cd =0;                %validation paper                     % non-dimensional damp
coef.Rho = 1000;              % fluid density
 
RadChan = 0.318;              %mono data
WidHull = 0.201;

 
 coef.A = WidHull*RadChan;          % wave channel area
 x0m = coef.GetDimCd('linear'); 
 x0p = -37.22/180*pi;     %validation paper             % non-dimensional damping coeff phase
 
 cdPlate = x0m*exp(1i*x0p);         % dimensional complex damping coefficient

% compute damped forces
FexD = waveChanDamping(Comp, CompP, cdPlate);
force_damped=abs(FexD);
 
%make power spectral density (index D for damping)           
Fs_D = flip(FexD.*A'); % force amplitude spectum
df = F0(2) - F0(1);
FD_D = 0.5*abs(Fs_D).^2./df; % force energy density spectrum
FexMod_D(n+j,:,:) = FD_D;
% surge force: DOF = 1
% heave force: DOF = 3
% pitch moment: DOF = 5

%plot all
plot(F0,FexMod(1,:,5),'b',F0,FexMod_D(1,:,5),'r','linewidth',1.5)
title('PSD(FexD) and PSD(FexU) for S2')

%plot(F0,force_damped(:,5),'b',F0,force_modelled(:,5),'r','linewidth',1.5)
%title('FexD and FexU for S2')


%plot exp and damped
%plot(F0,FexExp(2,:,3),'r',F0,FexMod_D(2,:,5),'b','linewidth',1.5)
%legend('Experimental', 'damped')

%plot difference between damped and undamped case
%plot (F0, (FexMod(1,:,5)-FexMod_D(1,:,5))/max(FexMod(1,:,5)))