% example script to apply the damping plate method for regular wave tests

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

%% load data
% direction = {'Fwd', 'Aft'}; % wave direction (0 pi)
% steepness = {'S120','S080'}; % wave steepness
% phi = [10 30 50]; % wave channel base plate angle

% example run
dataSetName = 'S120_Fwd_A30';

dataPath = 'C:\Users\Gabriel\OneDrive - Mocean Energy LTD\MACS-MIGS\ProjectData\experimental campaigns\2019\combined_monochromatic_data\'; % change this to your path
load([dataPath dataSetName '.mat']) % open the data

%% experimental data description

% inspect the experimental data

% force_experiment;  forces: dimension (wave frequency x DOF)
% phase_experiment; phases: dimension (wave frequency x DOF)

% surge force: DOF = 1
% heave force: DOF = 2
% pitch moment: DOF = 3

% T = opt.T; wave periods
% H = opt.H; water depth

%% simulation data (from WAMIT) description

% CompU;  undamped wave channel hydrodynamic computation
% Comp;   wave channel computation with damping plate
% Comp;   damping plate computation with wave channel

% NOTE: these are objects of FreqDomComp classes. Our mwave repository is required to work with these
% https://github.com/cmcnatt/mwave

%% inspect the undamped modelled data
Fex = squeeze(CompU.Fex);
force_modelled = abs(Fex); % magnitude of forces
phase_modelled = angle(Fex); % phase of forces from wave peak
        
 %% apply the damping plate method: try changing the coefficient
 
 % use ViscDampCoef class to dimensionalise old damping coeff
 coef = ViscDampCoef;
 
 coef.Cd = 2;                       % non-dimensional damping coeff magnitude
 coef.Rho = 1000;                   % fluid density
 coef.A = WidHull*RadChan;  % wave channel area
 x0m = coef.GetDimCd('linear'); 
 x0p = 70/180*pi;                   % non-dimensional damping coeff phase
 
 cdPlate = x0m*exp(1i*x0p);   % dimensional complex damping coefficient

 % compute damped forces
 FexD = waveChanDamping(Comp, CompP, cdPlate);
 
force_modelled_damped = abs(FexD); % magnitude of forces
phase_modelled_damped = angle(FexD); % phase of forces from wave peak

%% plot the pitch moment

%inds = [1 3 5]; % indices for surge, heave and pitch
figure;
subplot(1,2,1)
plot(1./T, force_experiment(:,3), 'ro', 1./T, force_modelled(:,5), 'k--', 1./T, force_modelled_damped(:,5), 'k')
xlabel('Frequency (Hz)')
ylabel( 'Pitch moment (Nm/m)')
legend('Experimental','Numerical undamped', 'Numerical damped')
subplot(1,2,2)
plot(1./T, phase_experiment(:,3), 'ro', 1./T, phase_modelled(:,5), 'k--', 1./T, phase_modelled_damped(:,5), 'k')
xlabel('Frequency (Hz)')
ylabel(' Phase (rad)')