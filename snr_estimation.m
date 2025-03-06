%% SNR Estimation
clc;
clear all;
close all;
%% Load data 
% Loads Reff, delay, wt_depth, beta, S, T 
load("..\PassiveSims\FinalSNRResults.mat")

%% Define Constants
Vsun = 5e-13;
Vsys = 9e-12;
df = 80e6; 
Tmax = 8;

%% Calculate SNR
% function of water table depth, sun elevation angle, and saturation/temperature
SNR = (Vsun^2/Vsys^2*df*Tmax).*R_eff;
SNR = (Vsun^2/Vsys^2*sqrt(df*Tmax)).*R_eff;

%% SNR Mean values, select scenarios for plotting
% % Varying saturation
% beta_1 = 30; % 25:45
% beta_kk = find(beta == beta_1); 
% 
% [X,Y] = meshgrid(S, wt_depth);
% Z = NaN(length(wt_depth), length(S));
% avgReff = NaN(length(wt_depth), length(S));
% for i = 1:length(wt_depth)
%     for j = 1:length(S)
%         Z(i,j) = nanmean(SNR(i,beta_kk,j,:),4);
%         avgReff(i,j) = mean(R_eff(i,beta_kk,j,:),4);
%     end
% end
% figure(1)
% pcolor(X,Y,10*log10(Z))
% shading interp
% c = colorbar;
% ylabel('Water table depth')
% xlabel('Saturation')
% ylabel(c,'SNR [dB]')
% title(["SNR Estimation", "Helheim Firn Aquifer", ...
%     sprintf("Elevation Angle %d deg", beta_1)])
% set ( gca, 'ydir', 'reverse' )
% 
% figure(2)
% pcolor(X,Y,avgReff)
% shading interp
% c = colorbar;
% ylabel('Water table depth')
% xlabel('Saturation')
% ylabel(c,'Average Reff')
% title(["Average Reff", "Helheim Firn Aquifer", ...
%     sprintf("Elevation Angle %d deg", beta_1)])
% 
% % Varying elevation
% T_1 = -0.001;
% T_kk = find(T == T_1,1)+1;
% 
% [X,Y] = meshgrid(beta, wt_depth);
% Z = NaN(length(wt_depth), length(beta));
% avgReff = NaN(length(wt_depth), length(beta));
% for i = 1:length(wt_depth)
%     for j = 1:length(beta)
%         Z(i,j) = nanmean(SNR(i,j, T_kk,:),4);
%         avgReff(i,j) = nanmean(R_eff(i,j,T_kk, :),4);
%     end
% end
% figure(3)
% pcolor(X,Y,10*log10(Z))
% shading interp
% c = colorbar;
% ylabel('Water table depth')
% xlabel('Elevation')
% ylabel(c,'SNR [dB]')
% title(["SNR Estimation", "Helheim Firn Aquifer", ...
%     sprintf("S/T %d / %d", T_1, S(T_kk))])
% set ( gca, 'ydir', 'reverse' )

%% SNR Median value, select scenarios for plotting
% Varying saturation
beta_1 = 30; % 25:45
beta_kk = find(beta == beta_1); 

[X,Y] = meshgrid(1:26, wt_depth);
Z = NaN(length(wt_depth), length(S));
medReff = NaN(length(wt_depth), length(S));
for i = 1:length(wt_depth)
    for j = 1:length(S)
        Z(i,j) = nanmedian(SNR(i,beta_kk,j,:),4);
        medReff(i,j) = nanmedian(R_eff(i,beta_kk,j,:),4);
    end
end
figure(4)
clf
hold on;
pcolor(X,Y,10*log10(Z))
contour(X, Y, Z, [10 10], 'r', 'LineWidth', 2);
shading interp
c = colorbar;
clim([0,max(10*log10(Z),[],'all')])
ylabel('Water table depth [m]')
TS = string(arrayfun(@(t, s) append(num2str(round(t,1)),' , ', ...
    num2str(round(s,2))), T, S, 'UniformOutput', false) );
xticks(1:2:26)
xticklabels(TS(1:2:26))
xlabel('Temperature [C], Saturation')
ylabel(c,'SNR [dB]')
title(["Median SNR Estimated", "Helheim Firn Aquifer", ...
    sprintf("Elevation Angle %d deg", beta_1)])
set ( gca, 'ydir', 'reverse' )

figure(5)
pcolor(X,Y,medReff)
shading interp
c = colorbar;
ylabel('Water table depth [m]')
TS = string(arrayfun(@(t, s) append(num2str(round(t,3)),' , ', ...
    num2str(round(s,2))), T, S, 'UniformOutput', false) );
xticks(1:2:26)
xticklabels(TS(1:2:26))
xlabel('Temperature [C], Saturation')
ylabel(c,'Median Reff')
title(["Median Reff", "Helheim Firn Aquifer", ...
    sprintf("Elevation Angle %d deg", beta_1)])

% Varying elevation
T_1 = -0.001;
T_kk = find(T == T_1,1)+1;

[X,Y] = meshgrid(beta, wt_depth);
Z = NaN(length(wt_depth), length(beta));
medReff = NaN(length(wt_depth), length(beta));
for i = 1:length(wt_depth)
    for j = 1:length(beta)
        Z(i,j) = nanmedian(SNR(i,j, T_kk,:),4);
        medReff(i,j) = nanmedian(R_eff(i,j,T_kk, :),4);
    end
end
figure(6)
clf
hold on;
pcolor(X,Y,10*log10(Z))
contour(X, Y, Z, [10 10], 'r', 'LineWidth', 2);
shading interp
c = colorbar;
ylabel('Water table depth [m]')
xlabel('Elevation [deg]')
ylabel(c,'SNR [dB]')
title(["Median SNR Estimated", "Helheim Firn Aquifer", ...
    sprintf("S/T %.2f / %.3f", round(S(T_kk),3), round(T_1,3))])
set ( gca, 'ydir', 'reverse' )

%% Finding Tmax
Vsun = 5e-13;
Vsys = 9e-12;
df = 80e6; 
SNR = 10; % 10dB = 10 linear

% Tmax = SNR/df*(Vsys^2/Vsun^2)./R_eff;
Tmax = 1/df*(SNR*(Vsys^2/Vsun^2)./R_eff).^2;


%% Select Scenarios for plotting
% Varying saturation
beta_1 = 30; % 25:45
beta_kk = find(beta == beta_1); 

[X,Y] = meshgrid(1:length(S), wt_depth);
Z = NaN(length(wt_depth), length(S));
for i = 1:length(wt_depth)
    for j = 1:length(S)
        Z(i,j) = nanmedian(Tmax(i,beta_kk,j,:),4);
    end
end
f = figure(10);
clf
hold on;
contour(X, Y, Z, [8 8], 'r', 'LineWidth', 2);
pcolor(X,Y,Z)
shading interp
c = colorbar;
ylabel('Water table depth [m]')
xlabel('Temperature [C], Saturation')
TS = string(arrayfun(@(t, s) append(num2str(round(t,1)),' , ', ...
    num2str(round(s,2))), T, S, 'UniformOutput', false) );
xticks(1:2:26)
xticklabels(TS(1:2:26))
ylabel(c,'Tmax [sec]')
set ( gca, 'ydir', 'reverse' )
title(["Tmax for SNR 10 dB", "Helheim Firn Aquifer", ...
    sprintf("Elevation Angle %d deg", beta_1)])

% Varying elevation
T_1 = -0.001;
T_kk = find(T == T_1,1)+1;

[X,Y] = meshgrid(beta, wt_depth);
Z = NaN(length(wt_depth), length(beta));
medReff = NaN(length(wt_depth), length(beta));
for i = 1:length(wt_depth)
    for j = 1:length(beta)
        Z(i,j) = nanmedian(Tmax(i,j, T_kk,:),4);
    end
end
figure(11)
clf
hold on;
pcolor(X,Y,Z)
contour(X, Y, Z, [8 8], 'r', 'LineWidth', 2);
shading interp
c = colorbar;
clim([0,1e0])
ylabel('Water table depth [m]')
xlabel('Elevation [deg]')
ylabel(c,'Tmax [sec]')
title(["Tmax for SNR 10 dB", "Helheim Firn Aquifer", ...
    sprintf("S/T %.2f / %.3f", round(S(T_kk),3), round(T_1,3))])
set ( gca, 'ydir', 'reverse' )