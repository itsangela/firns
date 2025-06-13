%% Multiple Heights Sun Bed Reflection Distance
%% Single Day
% Defining constants
hs = [5, 15, 25, 40]; % water table depth
freq = 330e6; % Center Frequency 330MHz
lambda = physconst('LightSpeed')/freq; % Wavelength
B = 15.36e6; % Bandwidth 15.36 MHz
n1 = 1; % refractive index of air

% refractive index of firn (n_ice = 1.78)
% for firn permitivity 2.49-5.49
n2 = sqrt(2.49); % sqrt(5.49)

location_name = 'Helheim Firn Aquifer, Greenland';

% Define time range and interval
tz = -2; % timezone
% Summer Solstice
% start_date = datetime([2026,06,21,0,0,0], ...
%     'Format','yyyy/MM/dd HH:mm:SS') - hours(tz);
% end_date = datetime([2026,06,22,0,1,0], ...
%     'Format','yyyy/MM/dd HH:mm:SS') -  hours(tz);
% % Spring Equinox
start_date = datetime([2026,03,20,0,0,0], ...
    'Format','yyyy/MM/dd HH:mm:SS') - hours(tz);
end_date = datetime([2026,03,20,23,59,0], ...
    'Format','yyyy/MM/dd HH:mm:SS') -  hours(tz);
time_interval = minutes(0:minutes(end_date-start_date));

% Define receiver position (lat, lon)
lat = 66.3575; % [66.353, 66.362] Bounding Coordinates
lon = -39.2235; % [-39.135, -39.312] Bounding Coordinates
altitude = 0;

% Convert datenum to UTC
mDateVec = start_date + time_interval;
UTC = string(mDateVec');

% Calculate AzEl
lat_vec = zeros(size(UTC,1),1)+lat;
lon_vec = zeros(size(UTC,1),1)+lon;
alt_vec = zeros(size(UTC,1),1)+altitude;
[sAz,sEl] = SolarAzEl(UTC, lat_vec, lon_vec, alt_vec);

belowHorizonkk = find(sEl<=0);
sEl(belowHorizonkk) = NaN;
sAz(belowHorizonkk) = NaN;

% Subsurface incidence angle
theta2 = asind(n1*sind(90-sEl)./n2);
% theta2 = asind(n1*cosd(sEl)./n2);

% Range
range = hs./cosd(theta2);

% horizontal distance
xs = 2*hs.*tand(theta2);
x_ref = xs./2.*cosd(sAz);
y_ref = xs./2.*sind(sAz);

figure(22)
clf
hold on;
p{1} = plot3(x_ref(:,1), y_ref(:,1), -ones(size(sAz))*hs(1), 'LineWidth',2);
p{2} = plot3(x_ref(:,2), y_ref(:,2), -ones(size(sAz))*hs(2), 'LineWidth',2);
p{3} = plot3(x_ref(:,3), y_ref(:,3), -ones(size(sAz))*hs(3), 'LineWidth',2);
p{4} = plot3(x_ref(:,4), y_ref(:,4), -ones(size(sAz))*hs(4), 'LineWidth',2);
start_kk = find(~isnan(y_ref(:,1)),1,'first');

% summer solstice (find index sun goes above horizon the second time)
% aboveH_kk = find(~isnan(y_ref(:,1)));
% start_kk = aboveH_kk(8);

plot3(x_ref(start_kk,1), y_ref(start_kk,1), -hs(1), 's', 'MarkerEdgeColor',"#0072BD", 'LineWidth',2);
plot3(x_ref(start_kk,2), y_ref(start_kk,2), -hs(2), 's', 'MarkerEdgeColor',"#D95319", 'LineWidth',2);
plot3(x_ref(start_kk,3), y_ref(start_kk,3), -hs(3), 's', 'MarkerEdgeColor',"#EDB120", 'LineWidth',2);
plot3(x_ref(start_kk,4), y_ref(start_kk,4), -hs(4), 's', 'MarkerEdgeColor',"#7E2F8E", 'LineWidth',2);
p{5} = plot3(nan,nan,nan,'ks');
p{6} = plot3(0, 0, 0, 'o','LineWidth',2, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
% summer solstice limits
% xlim([-25,40])
% % spring equinox limits
xlim([-40,25])
ylim([-40, 40])
xlabel('X distance (meters)')
ylabel('Y distance (meters)')
zlabel('Depth (meters)')
legend([p{:}], {'h=5m', 'h=15m', 'h=25m', 'h=40m', 'start', 'receiver'})
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d'));
title(['Map of Sun Bed Reflection Distance', location_name, date_string])
grid on;
