%% Maximum Integration Time

%% Single Day
% Defining constants
h = 40; % water table depth
freq = 330e6; % Center Frequency 330MHz
lambda = physconst('LightSpeed')/freq; % Wavelength
B = 15.36e6; % Bandwidth 15.36 MHz
n1 = 1; % refractive index of air
n2 = 1.5; % refractive index of firn (n_ice = 1.78)
location_name = 'Helheim Firn Aquifer, Greenland';

% Define time range and interval
tz = -2; % timezone
start_date = datetime([2025,06,21,0,0,0], ...
    'Format','yyyy/MM/dd HH:mm:SS') - hours(tz);
end_date = datetime([2025,06,22,00,0,0], ...
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
range = h./cosd(theta2);

% horizontal distance
xs = 2*h*tand(theta2);
x_ref = xs/2.*cosd(sAz);
y_ref = xs/2.*sind(sAz);

% Expected Delay
rd = 2*h.*tand(theta2).*cosd(sEl);
rice = 2*h./cosd(theta2); 
t_delay = rice./(physconst('LightSpeed')/n2) - rd./physconst('LightSpeed');

% Velocity
vel_ref = NaN(size(sAz));
for i =1:length(vel_ref)
    if i==1
        dx = x_ref(i+1)-x_ref(i);
        dy = y_ref(i+1)-y_ref(i);
        dr = sqrt(dx^2+dy^2);
        dt = minutes(mDateVec(i+1)-mDateVec(i));
        vel_ref(i) = dr/dt;
    elseif i < length(vel_ref)
        dx = x_ref(i+1)-x_ref(i-1);
        dy = y_ref(i+1)-y_ref(i-1);
        dr = sqrt(dx^2+dy^2);
        dt = minutes(mDateVec(i+1)-mDateVec(i-1));
        vel_ref(i) = dr/dt;
    else
        vel_ref(i) = vel_ref(i-1);
    end
end

vel_r = NaN(size(sAz));
for i =1:length(vel_r)
    if i < length(vel_r)
        % dtheta = theta2(i+1)-theta2(i);
        dtheta = deg2rad(sAz(i+1))-deg2rad(sAz(i));
        if dtheta < -2
            dtheta = dtheta+2*pi;
        end
        dt = minutes(mDateVec(i+1)-mDateVec(i));
        vel_r(i) = range(i)*dtheta/dt;
    else
        vel_r(i) = vel_r(i-1);
    end
end

% Fresnel zone
% r = sqrt(x_ref.^2+y_ref.^2)/2;
F = 2*sqrt(lambda*range);

% Integration time 
T_int = F./vel_ref;
T_int_r = F./vel_r;

% Ground Range Resolution
rg = (physconst('LightSpeed')/(2*n2*B))./sind(theta2);
% Ground Range Resolution vs sun angles 
sa = 0:41; % sun angles
rg_sa = (physconst('LightSpeed')/(2*n2*B))./(n1*sind(90-sa)./n2);

%Plot
figure(1)
clf
hold on;
plot(hours(time_interval), sEl, 'b', 'LineWidth',2);
plot(hours(time_interval), sAz, 'r', 'LineWidth',2);
% xticks(0:4:24)
xlim([0,24])
ylim([0, 360])
xlabel('Military Time Hours')
ylabel('Degrees')
legend('Elevation', 'Azimuth')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' local time';
title(['Sun Angles', location_name, date_string])
grid on

%Plot
figure(15)
clf
hold on;
plot(hours(time_interval), rad2deg(theta2), 'LineWidth',2);
% xticks(0:4:24)
xlim([0,24])
xlabel('Military Time Hours')
ylabel('Degrees')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' local time';
title(['Subsurface angle of incidence', location_name, date_string])
grid on

figure(2)
clf
subplot(1,2,1)
hold on;
plot(hours(time_interval), x_ref, 'b', 'LineWidth',2);
plot(hours(time_interval), y_ref, 'r', 'LineWidth',2);
plot(hours(time_interval), sqrt(x_ref.^2+y_ref.^2), 'k', 'LineWidth',2);
% xticks(0:4:24)
xlim([0,24])
legend('X position', 'Y position', '|R|')
xlabel('Military Time Hours')
ylabel('Meters')
title([sprintf('Map of Sun Bed Reflection Distances h=%d m', h), location_name, date_string])
subplot(1,2,2)
hold on;
plot(x_ref, y_ref, 'LineWidth',2);
% plot(x_ref(891), y_ref(891), 'x','LineWidth',2);
plot(0, 0, '*','LineWidth',2);
% xlim([min(x_ref),max(x_ref)])
% ylim([min(y_ref),max(y_ref)])
xlabel('X distance (meters)')
ylabel('Y distance (meters)')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' local time';
title([sprintf('Map of Sun Bed Reflection Distance h=%d m', h), location_name, date_string])
grid on

figure(3)
clf
subplot(3,1,1)
hold on;
plot(hours(time_interval), vel_ref, 'r','LineWidth',2);
% plot(hours(time_interval), vel_r, 'b--','LineWidth',2);
% xticks(0:4:24)
xlim([0,24])
xlabel('Military Time Hours')
ylabel('Meters/Minutes')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' local time';
title([sprintf('Velocity of Sun Azimuth Reflection Point h=%d m', h), location_name, date_string])
grid on
subplot(3,1,2)
hold on;
plot(hours(time_interval), F, 'b','LineWidth',2);
% xticks(0:4:24)
xlim([0,24])
xlabel('Military Time Hours')
ylabel('Meters')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' local time';
title([sprintf('Fresnel Zone h=%d m', h), location_name, date_string])
grid on
subplot(3,1,3)
hold on;
plot(hours(time_interval), T_int, 'k', 'LineWidth',2);
% xticks(0:4:24)
xlim([0,24])
xlabel('Military Time Hours')
ylabel('Minutes')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' local time';
title([sprintf('Integration Time h=%d m', h), location_name, date_string])
grid on

figure(4)
clf
subplot(1,2,1)
hold on;
plot(sa, rg_sa, 'LineWidth',2);
grid on;
xlim([0,41])
xlabel('Sun Angle (degrees)')
ylabel('Meters')
title(['Ground Range Resolution vs Sun Angle', location_name, date_string])
subplot(1,2,2)
hold on;
plot(hours(time_interval), rg, 'LineWidth',2);
xlabel('Military Time Hours')
ylabel('Meters')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' local time';
title(['Ground Range Resolution vs Time', location_name, date_string])
grid on

figure(5)
plot(hours(time_interval), t_delay.*1e9, 'LineWidth',2);
xlabel('Military Time Hours')
ylabel('Delay (nanosec)')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' local time';
title([sprintf('Expected Delay Time h=%d m', h), location_name, date_string])
grid on

%% Multiple Heights
% Defining constants
freq = 330e6; % Center Frequency 330MHz
lambda = physconst('LightSpeed')/freq; % Wavelength
n1 = 1; % refractive index of air
n2 = 1.5; % refractive index of firn (n_ice = 1.78)

% Define time range and interval
tz = -2; % timezone
start_date = datetime([2025,03,01,0,0,0], ...
    'Format','yyyy/MM/dd HH:mm:SS') - hours(tz);
end_date = datetime([2025,03,01,23,59,0], ...
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

belowHorizonkk = find(sEl<0);
sEl(belowHorizonkk) = NaN;
sAz(belowHorizonkk) = NaN;

% Subsurface incidence angle
theta2 = asin(n1*sind(90-sEl)./n2);

% Defining heights
water_table_heights = [10,20,30,40];
figure(6)
clf
for h = water_table_heights
    % Range
    range = h./cos(theta2);
    
    % horizontal distance
    xs = 2*h*tan(theta2);
    x_ref = xs/2.*cosd(sAz);
    y_ref = xs/2.*sind(sAz);
    
    % Velocity
    vel_ref = NaN(size(sAz));
    for i =1:length(vel_ref)
        if i < length(vel_ref)
            dx = x_ref(i+1)-x_ref(i);
            dy = y_ref(i+1)-y_ref(i);
            dr = sqrt(dx^2+dy^2);
            dt = minutes(mDateVec(i+1)-mDateVec(i));
            vel_ref(i) = dr/dt;
        else
            vel_ref(i) = vel_ref(i-1);
        end
    end
    
    % Fresnel zone
    F = 2*sqrt(lambda*range);
    
    % Integration time 
    T_int = F./vel_ref;
    
    %Plot
    % figure(4)
    % clf
    % hold on;
    % plot(hours(time_interval), sEl, 'b', 'LineWidth',2);
    % plot(hours(time_interval), sAz, 'r', 'LineWidth',2);
    % xticks(0:4:24)
    % xlim([0,24])
    % xlabel('Military Time Hours')
    % ylabel('Degrees')
    % legend('Elevation', 'Azimuth')
    % date_string = string(datetime(start_date, 'Format', 'yyyy MMM d'));
    % title(['Sun Angles', location_name, date_string])
    % grid on
    
    % figure(5)
    % clf
    % subplot(1,2,1)
    % hold on;
    % plot(hours(time_interval), x_ref, 'b', 'LineWidth',2);
    % plot(hours(time_interval), y_ref, 'r', 'LineWidth',2);
    % plot(hours(time_interval), sqrt(x_ref.^2+y_ref.^2), 'k', 'LineWidth',2);
    % xticks(0:4:24)
    % xlim([0,24])
    % legend('X position', 'Y position', '|R|')
    % xlabel('Military Time Hours')
    % ylabel('Meters')
    % title(['Map of Sun Bed Reflection Distances h=40m', location_name, date_string])
    % subplot(1,2,2)
    % hold on;
    % plot(x_ref, y_ref, 'LineWidth',2);
    % plot(0, 0, '*','LineWidth',2);
    % xlabel('X distance (meters)')
    % ylabel('Y distance (meters)')
    % date_string = string(datetime(start_date, 'Format', 'yyyy MMM d'));
    % title(['Map of Sun Bed Reflection Distance h=40m', location_name, date_string])
    % grid on
    
    subplot(3,1,1)
    hold on;
    plot(hours(time_interval), vel_ref ,'LineWidth',2);
    xticks(0:4:24)
    xlim([0,24])
    xlabel('Military Time Hours')
    ylabel('Meters/Minutes')
    date_string = string(datetime(start_date, 'Format', 'yyyy MMM d'));
    title(['Velocity of Sun Azimuth Reflection Point', location_name, date_string])
    grid on
    subplot(3,1,2)
    hold on;
    plot(hours(time_interval), F ,'LineWidth',2);
    xticks(0:4:24)
    xlim([0,24])
    xlabel('Military Time Hours')
    ylabel('Meters')
    date_string = string(datetime(start_date, 'Format', 'yyyy MMM d'));
    title(['Fresnel Zone', location_name, date_string])
    grid on
    subplot(3,1,3)
    hold on;
    plot(hours(time_interval), T_int, 'LineWidth',2);
    xticks(0:4:24)
    xlim([0,24])
    xlabel('Military Time Hours')
    ylabel('Minutes')
    date_string = string(datetime(start_date, 'Format', 'yyyy MMM d'));
    title(['Integration Time', location_name, date_string])
    grid on
end
subplot(3,1,1)
legend('h=10 m', 'h=20 m', 'h=30 m', 'h=40 m')


%% Multiple Permitivities
% Defining constants
h = 40;
freq = 330e6; % Center Frequency 330MHz
lambda = physconst('LightSpeed')/freq; % Wavelength
n1 = 1; % refractive index of air

% Define time range and interval
tz = -2; % timezone
start_date = datetime([2025,03,01,0,0,0], ...
    'Format','yyyy/MM/dd HH:mm:SS') - hours(tz);
end_date = datetime([2025,03,01,23,59,0], ...
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

belowHorizonkk = find(sEl<0);
sEl(belowHorizonkk) = NaN;
sAz(belowHorizonkk) = NaN;

% Defining heights
n2_s = linspace(sqrt(2.49), sqrt(5.9), 5);
figure(6)
clf
for n2 = n2_s
    % Subsurface incidence angle
    theta2 = asin(n1*sind(90-sEl)./n2);
    
    % Range
    range = h./cos(theta2);
    
    % horizontal distance
    xs = 2*h*tan(theta2);
    x_ref = xs/2.*cosd(sAz);
    y_ref = xs/2.*sind(sAz);
    
    % Velocity
    vel_ref = NaN(size(sAz));
    for i =1:length(vel_ref)
        if i < length(vel_ref)
            dx = x_ref(i+1)-x_ref(i);
            dy = y_ref(i+1)-y_ref(i);
            dr = sqrt(dx^2+dy^2);
            dt = minutes(mDateVec(i+1)-mDateVec(i));
            vel_ref(i) = dr/dt;
        else
            vel_ref(i) = vel_ref(i-1);
        end
    end
    
    % Fresnel zone
    F = 2*sqrt(lambda*range);
    
    % Integration time 
    T_int = F./vel_ref;
    
    subplot(3,1,1)
    hold on;
    plot(hours(time_interval), vel_ref ,'LineWidth',2);
    xticks(0:4:24)
    xlim([0,24])
    xlabel('Military Time Hours')
    ylabel('Meters/Minutes')
    date_string = string(datetime(start_date, 'Format', 'yyyy MMM d'));
    title([sprintf('Velocity of Sun Azimuth Reflection Point h=%d m', h), ...
        location_name, date_string])
    grid on
    subplot(3,1,2)
    hold on;
    plot(hours(time_interval), F ,'LineWidth',2);
    xticks(0:4:24)
    xlim([0,24])
    xlabel('Military Time Hours')
    ylabel('Meters')
    date_string = string(datetime(start_date, 'Format', 'yyyy MMM d'));
    title([sprintf('Fresnel Zone h=%d m', h), location_name, date_string])
    grid on
    subplot(3,1,3)
    hold on;
    plot(hours(time_interval), T_int, 'LineWidth',2);
    xticks(0:4:24)
    xlim([0,24])
    xlabel('Military Time Hours')
    ylabel('Minutes')
    date_string = string(datetime(start_date, 'Format', 'yyyy MMM d'));
    title([sprintf('Integration Time h=%d m', h), location_name, date_string])
    grid on
end
subplot(3,1,1)
legend('n2 = 1.5780', 'n2 = 1.7907', 'n2 = 2.0035','n2 = 2.2162', 'n2=2.4290')