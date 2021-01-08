%% load and plot monochromatic wave data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc


direction = {'Fwd', 'Aft'};     % wave direction [0, pi]
steepness = {'S120','S080'};    % wave steepness height/wavelength [1/120, 1/80]
phi = [10 30 50];               % wave channel base plate angle

for j =1:2 % loop over each wave direction
    
    for k =1:2 % loop over each wave steepness
        
        
        for n = 1:length(phi)
            clear('opt')
            dataPath = 'C:\Users\Gabriel\OneDrive - Mocean Energy LTD\MACS-MIGS\ProjectData\experimental campaigns\2019\combined_monochromatic_data\'; % change this to your path
            dataSetName = [steepness{k} '_' direction{j} '_A' num2str(phi(n))]; % set containing experimental and numerical data
            load([dataPath dataSetName '.mat']) % open the data
            
            % inspect the experimental data
            force_experiment = opt.ExForces; % forces: dimension (wave frequency x DOF)
            phase_experiment = opt.ExPhases; % phases: dimension (wave frequency x DOF)
            % surge force: DOF = 1
            % heave force: DOF = 2
            % pitch moment: DOF = 3
            T = opt.T; % wave periods
            H = opt.H; % water depth
            
            % inspect the simulation data (from WAMIT)
            CompU = opt.CompU{1};   % undamped wave channel hydrodynamic computation
            Comp = opt.Comp{1};     % wave channel computation with damping plate
            CompP = opt.CompP{1};   % damping plate computation with wave channel
            
            % NOTE: these are objects of (FreqDomComp) classes. Our mwave repository is required to work with these
            % https://github.com/cmcnatt/mwave
            
            % inspect the undamped modelled forces
            
            Fex = squeeze(CompU.Fex);
            force_modelled = abs(Fex);   % magnitude of forces
            phase_modelled = angle(Fex); % phase of forces from wave peak
            % surge force: DOF = 1
            % heave force: DOF = 3
            % pitch moment: DOF = 5
            
            % plot the data
            ylab = {'Surge force (N/m)', 'Heave force (N/m)', 'Pitch moment(Nm/m)'};
            ylab2 = {'Phase (rad)', 'Phase (rad)', ' Phase (rad)'};
            inds = [1 3 5]; % indices for surge, heave and pitch
            figure;
            a=0;
            for m = 1:3
                subplot(3,2,m+a)
                plot(1./T, force_experiment(:,1), 'ro', 1./T, force_modelled(:,inds(m)), 'k--')
                xlabel('Frequency (Hz)')
                ylabel(ylab{m})
                legend('Experimental','Numerical')
                subplot(3,2,2*m)
                plot(1./T, phase_experiment(:,1), 'ro', 1./T, phase_modelled(:,inds(m)), 'k--')
                xlabel('Frequency (Hz)')
                ylabel(ylab2{m})
                legend('Experimental','Numerical')
                a=a+1;
            end
        end % phi
        
    end % steepness
end % direction

