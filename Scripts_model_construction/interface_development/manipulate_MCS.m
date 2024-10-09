clear all, close all, clc

out_dir   = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/';


theMCS    = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/interface_depthConv_manipulated.txt';
MCS_oc_tb = table2array(readtable(theMCS));

% imported MCS data
MCS_lon   = MCS_oc_tb(:,1);
MCS_lat   = MCS_oc_tb(:,2);
MCS_depth = MCS_oc_tb(:,3);


figure(98453), clf
plot3(MCS_lon, MCS_lat, MCS_depth, '.b')
grid on



% imported inverted points
% % add inverted points -- TEMP SOLUTION
invP_tb   = table2array(readtable('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/interface_invertedPoints.txt')); 
invP_lon  = invP_tb(:,1);
invP_lat  = invP_tb(:,2);
invP_dp   = invP_tb(:,3)+6;

MCS_inv_lon   = vertcat(MCS_lon, invP_lon);
MCS_inv_lat   = vertcat(MCS_lat, invP_lat);
MCS_inv_depth = vertcat(MCS_depth, invP_dp);




figure(45), clf
plot(MCS_inv_lon, MCS_inv_lat, '.r')
grid on
title('select a subset of shots with four selected points')
[x,y] = ginput(4);

lat_ind       = find(max(y)>MCS_inv_lat & min(y)<MCS_inv_lat);
mcs_lon_filt1 = MCS_inv_lon(lat_ind); mcs_lat_filt1 = MCS_inv_lat(lat_ind); mcs_dp_filt1 = MCS_inv_depth(lat_ind);
lon_ind       = find(max(x)>mcs_lon_filt1 & min(x)<mcs_lon_filt1);
mcs_lon_filt2 = mcs_lon_filt1(lon_ind); mcs_lat_filt2 = mcs_lat_filt1(lon_ind); mcs_dp_filt2 = mcs_dp_filt1(lon_ind);

figure(46),clf
plot(MCS_inv_lon, MCS_inv_lat, '.r')
hold on
plot(mcs_lon_filt2, mcs_lat_filt2, '.b')
grid on
hold on
title('First Two Points')
[x1, y1] = ginput(2);
plot(x1, y1, '-m')
hold on
title('Second Two Points')
[x2, y2] = ginput(2);
plot(x2, y2, '-m')

% first extended line
x1_ext = linspace(min(x1), max(x1), 25);
y1_ext = linspace(min(y1), max(y1), 25);
z1_ext = zeros(length(x1_ext))+mean(mcs_dp_filt2); z1_ext = z1_ext(:,1);


% second extended line
x2_ext = linspace(min(x2), max(x2), 25);
y2_ext = linspace(min(y2), max(y2), 25);
z2_ext = zeros(length(x2_ext))+mean(mcs_dp_filt2); z2_ext = z2_ext(:,1);


mcs_lon_combo   = vertcat(MCS_inv_lon,   x1_ext', x2_ext');
mcs_lat_combo   = vertcat(MCS_inv_lat,   y1_ext', y2_ext');
mcs_depth_combo = vertcat(MCS_inv_depth, z1_ext, z2_ext);

figure(47), clf
plot3(mcs_lon_combo, mcs_lat_combo, mcs_depth_combo, '.r')
grid on


%% Write

fid = fopen(append(append(out_dir, 'interface_depthConv_manipulated.txt')), 'w');

if fid == -1
    error('Unable to opne file.')
end

for i = 1:length(mcs_lon_combo)
    fprintf(fid, '%f %f %f\n', mcs_lon_combo(i), mcs_lat_combo(i), mcs_depth_combo(i));
end

fclose(fid);
