% WGN_Firn

% Use Parameters 

% clear all; close all; clc; 

%======================   CONTROL PANEL   ============================
% We can control the chosen parameter sweep by the following numbers
%
addpath('TimeSeries')
%load TimeSeriesOutput; % load('Scenario2'); load('Scenario1');
load TimeSeries_Wideband

% TimeSeries_Wideband.mat
	% > delay - 1x562 vector - delay between direct path and reflection from the water table in seconds
	% > density - 2497x1 vector - FA16-4 density profile as a function of depth in kg/m^3, average density profile from the SUMup2022 collection, ice 
	%   lenses have been added back in based on the stratigraphy published in Miller et al. (2020)
	% > depth - 2497x1 vector - depth below the surface in meters
	% > fc - 1x161 vector - frequency steps across the whole bandwidth in Hz used for the reflectivity and loss simulations
	% > loss - 161x562 vector (frequency x time) - total power lost due to transmission and absorption in the overlying firn in linear power units as 
	%   a function of frequency and time
	% > R - 161x562 vector (frequency x time) - power reflection coefficient of the water table in linear power units as a function of frequency and time
	% > R_eff - 161x562 vector (frequency x time)- effective power reflection coefficient (combines R and loss) in linear power units as a function of 
	%   frequency and time
	% > radar_depth - 1x562 vector - depth of the radar system below the surface in meters, based on RACMOv2.3p2 estimates of daily snow accumulation
	% > S - 2497x562 matrix - saturation of the subsurface as a function of depth and time
	% > sEl - 1x562 vector - sun elevation angle in degrees as a function of time
	% > temperature - 2497x562 matrix - firn temperature  in degrees C as a function of depth and time, interpolated from the FA15-1 thermistor 
	%   measurements with corrections for change in surface elevation based on RACMOv2.3p2 estimates of daily snow accumulation
	% > UTC - 1x562 vector - date in MATLAB datetime format
	% > wt_depth - 1x562 vector - water table depth below the surface in meters as a function of time, corrected for change in surface elevation 
	%   estimated from RACMOv2.3p2 daily snow accumulation, not accounting for firn compaction

no1 = 1;     % number of Simulations -  % 65000; % 1000  % small number to get idea of what it looks like

% ==================  Graphing CONTROL PANEL   =======================

plotLAGS = 0; % if you want to do the plot ''plotLAGS'' then set to 1, otherwise 0;

%======================   CONTROL PANEL   ============================

set(0,'defaultLineLineWidth',2);   set(0,'defaultAxesFontSize', 18);

c = 2.998*10^8;                      % m/s speed of light
f = 330e6;                           % center frequency
%T = 100e-3;%0.01;%100e-6;                    % integration (receiving) time in s

T = 0.1;
%T = 10000e-6;                    % integration (receiving) time in s
Fs = 80e6;                           % 15.36e6;     % sampling rate in Hz
%Fs = 15.36e6;

N = T*Fs                            % total number of samples received
Nmin = min(10e3,N);                  % minimum samples needed to see example echo
fVals = linspace(-Fs/2,Fs/2, N).';   % frequency domain spacing

% R_eff - 161x562 vector (frequency x time)- effective power reflection coefficient 
% (combines R and loss) in linear power units as a function of frequency and time
EffReflect = R_eff; 

% 1x562 vector - delay between direct path and water table reflection  in sec
DelayTime = delay; 

jj = 1;
N1 = length(EffReflect(:,jj));
fVals1 = linspace(-Fs/2,Fs/2, N1).';   % frequency domain spacing correction


%
nSims = size(EffReflect,2); % no1; % 65000; % small number to get idea of what it looks like

testindx = 1;
%alpha0 = ElevationAngle(testindx)*180/pi;
%for ii = 430:480
%for ii = 340:490

xSim = 1:nSims;
for ii = xSim
    
    Rs = EffReflect(:,ii);   % individual reflectivity band
    Rs_interp = interp1(fVals1,Rs,fVals); % 

    tau_refl = DelayTime(ii);
    
    % generate the white noise
    y = wgn(N,1,(4*10^(-14))*(Fs)/50,50,'linear','complex');      % Sun
    n = wgn(N,1,(1.4*10^(-13))*(Fs)/50,50,'linear','complex');    % Noise
    
    % Construct White Noise Echo Signal
    F_wgn = fft(y);
    F_d = F_wgn;                                        % DIRECT
    F_r = Rs_interp.*F_d.*exp(1j*(-2*pi*(tau_refl)*fVals));    % REFLECTED 
    F_bn = fft(n);                                      % BACKGROUND NOISE

    F = F_d + F_r + F_bn; % direct + reflected + noise

        % Autocorrelation via FFT computations (Wiener-Khinchin Theorem)
    PSD = (1/(Fs*N))*abs(F).^2;   % PSD
    Rxx = ifft(PSD);              % Autocorrelation
    
    Rxx = Rxx/max(abs(Rxx));    % normalize
    
    Rxx_p = Rxx(1:Nmin);                 % truncate
    samp_num = round(tau_refl*80e6+1);

    peak = Rxx_p(samp_num);
    
    lags = (0:length(Rxx_p)-1)/Fs;       % delay times

    Rnorm = abs(Rxx_p)/max(abs(Rxx_p));  % absolute magnitude
    % Rnorm = abs(Rxx_p)/max(abs(Rxx_p(2:1536)));  % absolute magnitude
    Rnorm_sq = Rnorm.^2;                 % square for power
    
    
    %xmin = max(tau_refl)*2.5*1e6;           % xlim in microsec
    xmin = 1.2;

    if(plotLAGS)
        figure(1), set(gcf,'color','w');
        %hold on
        plot(lags*1e6, 10*log10(Rnorm_sq),'k')
        xlim([0 xmin]), ylim([-60 0])
        ylabel('Magnitude (dB)')
        xlabel('Delay (\musec)')
        %title(['Autocorrelation, t:',num2str(T*1e3),'ms, d: ', num2str(d), 'm, \alpha: ',num2str(alpha0),'\circ'])
        title('Autocorrelation')
        hold on
        xline(tau_refl*1e6,'r')
        hold off
        pause(0.01)
    end
    
    % F_save(:,jj) = F;

    Rxx_peak(jj) = peak;

    Rxx_phase(jj) = phase(peak);
    Rxx_angle(jj) = angle(peak);

    Rnorm_save(:,jj) = Rnorm;
    Rnorm_sq_save(:,jj) = Rnorm_sq;
    Rnorm_sq_save_dB(:,jj) = 10*log10(Rnorm_sq);
    jj = jj+1;
ii

end

saveData = 0;
if(saveData)
    save("PassiveFirn_N8000_FS10MHz.mat","F_save","fVals","T","Fs","N", "Nmin","Rnorm_save","lags")
    save("PassiveFirnParameters.mat","fVals","T","Fs","N", "Nmin","Rnorm_save","lags")
    %writematrix(F_save,'PassiveFirnFreqDomain.csv')
end

%% Plotting echo power with delay
figure(1)
clf
[X,Y] = meshgrid(1:562, lags*1e6);
Z = Rnorm_sq_save_dB;
pcolor(X,Y,Z)
xlim([1,562])

dates_to_label = find(day(UTC) ==1);
xticks(dates_to_label)
xticklabels(string(datetime(UTC(dates_to_label),'Format', 'd MMM yyyy')))
xtickangle(45)

shading interp
cb = colorbar;
clim([-60,0])
ylim([0,0.45])
ylabel('Delay [microsec]')
xlabel('Date')
ylabel(cb,'Normalized Power [dBW]')
title(["Time Series", "10 ms Integration Time, 80 MHz BW"])
set ( gca, 'ydir', 'reverse' )
set(gca,'Color','black')
set(gca,'XMinorTick','off','YMinorTick','off')
set(gca,'TickDir','out')
box off

%% Plotting echo power with height
figure(2)
clf
[X,Y] = meshgrid(1:562, lags*c);
Z = Rnorm_sq_save_dB;
pcolor(X,Y,Z)
xlim([1,562])

dates_to_label = find(day(UTC) ==1);
xticks(dates_to_label)
xticklabels(string(datetime(UTC(dates_to_label),'Format', 'd MMM yyyy')))
xtickangle(45)

shading interp
cb = colorbar;
clim([-60,0])
ylim([0,0.45e-6*c])
ylabel('Height [m]')
xlabel('Date')
ylabel(cb, 'Normalized Power [dBW]')
title(["Time Series", "10 ms Integration Time, 80 MHz BW"])
set ( gca, 'ydir', 'reverse' )
set(gca,'Color','black')
set(gca,'XMinorTick','off','YMinorTick','off')
set(gca,'TickDir','out')
box off