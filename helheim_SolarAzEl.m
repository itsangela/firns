%% Solar Az El for Helheim Firm Aquifer

%% Single Day
% Define time range and interval
start_date = datetime([2025,03,01,0,0,0],'Format','yyyy/MM/dd HH:mm:SS');
end_date = datetime([2025,03,01,23,59,0],'Format','yyyy/MM/dd HH:mm:SS');
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

%Plot
figure(1)
clf
hold on;
plot(hours(time_interval), sEl, 'b', 'LineWidth',2);
plot(hours(time_interval), sAz, 'r', 'LineWidth',2);
xticks(0:4:24)
xlim([0,24])
xlabel('Military Time Hours')
ylabel('Degrees')
legend('Elevation', 'Azimuth')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d'));
title(['Sun Angles', 'Helheim Firn Aquifer', date_string])
grid on

%% Month
% Define time range and interval
start_date = datetime([2025,03,01,0,0,0],'Format','yyyy/MM/dd HH:mm:SS');
end_date = datetime([2025,03,31,23,0,0],'Format','yyyy/MM/dd HH:mm:SS');
time_interval = hours(0:hours(end_date-start_date));

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

% Plot
figure(2)
clf
hold on;
plot(mDateVec, sEl, 'b');
plot(mDateVec, sAz, 'r');
xlim([mDateVec(1),mDateVec(end)])
xlabel('Dates')
ylabel('Degrees')
legend('Elevation', 'Azimuth')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' - ' + string(datetime(end_date, 'Format', 'yyyy MMM d'));
title(['Sun Angles', 'Helheim Firn Aquifer', date_string])
grid on;

%% Summer Season
% Define time range and interval
start_date = datetime([2025,03,01,0,0,0],'Format','yyyy/MM/dd HH:mm:SS');
end_date = datetime([2025,09,30,23,0,0],'Format','yyyy/MM/dd HH:mm:SS');

% 24 hrs every 15 days
time_interval = [];
plt_date_vec = []; % for plot
plt_date_labels = [];
for dd = days(0:15:days(end_date-start_date))
    time_interval = [time_interval, dd + hours(0:23)];
    plt_date_vec = [plt_date_vec, days(dd/15+ hours(0:23))];
    plt_date_labels = [plt_date_labels, start_date+dd];
end

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

% Subsurface incidence angle
theta2 = asin(n1*sind(90-sEl)./n2);

% Plot
figure(3)
clf
hold on;
plot(plt_date_vec, sEl, 'b');
plot(plt_date_vec(1:24:end), sEl(1:24:end), 'b*');
plot(plt_date_vec, sAz, 'r');
plot(plt_date_vec(1:24:end), sAz(1:24:end), 'r*');
xticks(0:floor(days(end_date-start_date)))
xticklabels(string(datetime(plt_date_labels,'Format', 'MMM d')))
% xlim([0,24])
xlabel('Dates')
ylabel('Degrees')
legend('Elevation', 'Azimuth')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' - ' + string(datetime(end_date, 'Format', 'yyyy MMM d'));
title(['Sun Angles', 'Helheim Firn Aquifer', date_string])
grid on

%%
figure(4)
clf
hold on;
plot(plt_date_vec, rad2deg(theta2), 'b');
plot(plt_date_vec(1:24:end), rad2deg(theta2(1:24:end)), 'b*');
xticks(0:floor(days(end_date-start_date)))
xticklabels(string(datetime(plt_date_labels,'Format', 'MMM d')))
% xlim([0,24])
xlabel('Dates')
ylabel('Degrees')
date_string = string(datetime(start_date, 'Format', 'yyyy MMM d')) + ...
    ' - ' + string(datetime(end_date, 'Format', 'yyyy MMM d'));
title(['Subsurface Angle of Incidence', 'Helheim Firn Aquifer', date_string])
grid on