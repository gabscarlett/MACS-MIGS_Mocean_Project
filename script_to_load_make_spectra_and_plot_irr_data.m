% script to make spectra from iregular wave experimental time series
% and reqular modelled frequency domain data

% for the irregular wave tests we need to make measured and modelled force spectra to compare

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

for k = 1 % forward direction
    
    for m =1:length(Tp)
        
        for n=1:length(phi)
            
            
            %% load experimental data
            dataPath = 'C:\Users\Gabriel\OneDrive - Mocean Energy LTD\MACS-MIGS\ProjectData\experimental campaigns\2019\UoE Flume Tests 2019-02\IrrTS';
            modelledDataPath = 'C:\Users\Gabriel\OneDrive - Mocean Energy LTD\MACS-MIGS\ProjectData\experimental campaigns\2019\modelled_irregular_data\';
            dataSetName = [direction{k} '_A' num2str(phi(n)) '_Irr00TS' num2str(m)];
            forceTab = readtable([dataPath '\' dataSetName '.csv']);
            
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
            
            % frequency from 0.5 Hz to 1.5 Hz vector to line up with freqs
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
            
            %% get the modelled spectra
            
            % load modelled data (run over F0)
            
            load([modelledDataPath dataSetName '_modelled']);
            
            % get the undamped forces
            FexU = squeeze(CompU.Fex);            
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
        end
        j=j+3;
        
    end
    
end


%% plot the pitch moment spectra for each data set

% ****NOTE for modelled data:  inds = [1 3 5]; for surge, heave and pitch

figure;
for n = 1:18
    subplot(6,3,n);
    plot(F0, FexExp(n,:,3),'r',F0, FexMod(n,:,5),'k--');
    xlim([0 2]);
    xlabel('f (Hz)') 
    ylabel('Pitch (N^2 m^2 Hz^{-1})')
    ntitle(myCell(n),'location','northeast','fontsize',9)
        if n==1      
        legend({'data', 'Modelled'}.....
        ,'Location','northwest','NumColumns',2)
    end
end