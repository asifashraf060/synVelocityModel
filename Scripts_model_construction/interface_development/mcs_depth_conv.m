%% SCRIPT: CONVERSION OF MCS Travel Time to Depth using inverted Velocity model
% Asif, Jul 14, 2023
clear all, close all, clc
%% INPUT
theMCS    = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/time_TOC_SFsubtracted.txt';
theModel  = '/Users/asifashraf/Talapas/tlOutput/231025_112721/srModel_it7.mat';   
out_dir   = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/';

%% CONSISTENT INPUT

%column for longitude values in the table
lon_col   = 1;
%column for latitude values in the table
lat_col   = 2;
%column for depth values of the slab
depth_col = 3;

%% load input files

load(theModel)

MCS_oc_tb = table2array(readtable(theMCS));

% imported MCS data
MCS_lon   = MCS_oc_tb(:,lon_col);
MCS_lat   = MCS_oc_tb(:,lat_col);
MCS_time  = (MCS_oc_tb(:,depth_col));

figure(11), clf
plot3(MCS_lon, MCS_lat, MCS_time, '.b')

avg_vs = [];

%% Calculation

%% Convert
avg_vs = [];
evs    = [];
for i = 1:length(MCS_lon)
   
   % convert lat, lon to experiment x, y
   [x, y] = map2xy(MCS_lon(i), MCS_lat(i), srModel.srGeometry);
   % find the index for that x and y
   ind_x  = find(abs(srModel.xg-x) == min(abs(srModel.xg-x)));
   ind_y  = find(abs(srModel.yg-y) == min(abs(srModel.yg-y)));
   % extract the velocity in that x,y point
   u_1D   = squeeze(srModel.P.u(ind_x ,ind_y ,:));
   
   % get the elevation at that xy point
   ev     = srModel.elevation(ind_x, ind_y);

   t_1D   = abs(srModel.zg).*u_1D;
   t_sum  = cumsum(t_1D);
   
   

   % find where that depth value is
   ind_z  = find(abs(t_sum-MCS_time(i)) == min(abs(t_sum-MCS_time(i))));
   % average of the velocities up to that depth
   avg_u  = mean(u_1D(1:ind_z));
   avg_v  = 1/avg_u;

   avg_vs(i,:) = avg_v;
   evs(i,:)    = ev;
end

avg_vs = imgaussfilt(avg_vs, 10);


MCS_depth = (MCS_time.*avg_vs) + evs;

MCS_depth_filt = imgaussfilt(MCS_depth, 10);


figure(1),clf
plot3(MCS_lon, MCS_lat, (MCS_depth), '.r')
grid on
%set(gca, 'zDir', 'reverse')
hold on
plot3(MCS_lon, MCS_lat, (MCS_depth_filt), '.b')
plot3(srModel.LON, srModel.LAT, srModel.elevation);
scatter(MCS_lon, MCS_lat, 25 , avg_vs, 'filled')
cc = colorbar;
cc.Label.String = 'Average Vp (km/s)';
xlabel('Longitude'), ylabel('Latitude'), zlabel('Depth (km)')
set(gca, 'FontSize', 16)

%% Write

fid = fopen(append(append(out_dir, 'interface_depthConv.txt')), 'w');

if fid == -1
    error('Unable to opne file.')
end

for i = 1:length(MCS_time)
    fprintf(fid, '%f %f %f %f\n', MCS_lon(i), MCS_lat(i), MCS_time(i), MCS_depth(i));
end

fclose(fid);


