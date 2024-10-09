clear all, close all, clc

out_dir   = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/';


theMCS    = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/interface_depthConv_manipulated.txt';
MCS_oc_tb = table2array(readtable(theMCS));

%% IMPORT

% imported MCS data
MCS_lon   = MCS_oc_tb(:,1);
MCS_lat   = MCS_oc_tb(:,2);
MCS_depth = MCS_oc_tb(:,3);


% imported inverted points
% % add inverted points -- TEMP SOLUTION
% invP_tb   = table2array(readtable('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/interface_invertedPoints.txt')); 
% invP_lon  = invP_tb(:,1);
% invP_lat  = invP_tb(:,2);
% invP_dp   = invP_tb(:,3)+6;

MCS_inv_lon   = vertcat(MCS_lon); %invP_lon);
MCS_inv_lat   = vertcat(MCS_lat); %invP_lat);
MCS_inv_depth = vertcat(MCS_depth); %invP_dp);


%% SELECTION

figure(45), clf
plot(MCS_inv_lon, MCS_inv_lat, '.r')
grid on
title('select a subset of shots with four selected points')
[x,y] = ginput(4);

lat_ind       = find(max(y)>MCS_inv_lat & min(y)<MCS_inv_lat);
lon_ind       = find(max(x)>MCS_inv_lon & min(x)<MCS_inv_lon);
idx           = intersect(lat_ind, lon_ind);

figure(45), clf
plot(MCS_inv_lon, MCS_inv_lat, '.r')
hold on
plot(MCS_inv_lon(idx), MCS_inv_lat(idx), '.b')
title('Check the selected shots in blue || Enter if everything is ok')
ij = ginput();



%% DELETION
if isempty(ij)
    MCS_inv_lon(idx) = [];
    MCS_inv_lat(idx) = [];
    MCS_inv_depth(idx) = [];
    
    figure(45), clf
    plot(MCS_inv_lon, MCS_inv_lat, '.r')
    hold on
    title('Check the deletion')
end


figure(46),clf
plot3(MCS_inv_lon, MCS_inv_lat, MCS_inv_depth, '.r')
hold on
grid on



%% Write

fid = fopen(append(append(out_dir, 'interface_depthConv_manipulated.txt')), 'w');

if fid == -1
    error('Unable to opne file.')
end

for i = 1:length(MCS_inv_lon)
    fprintf(fid, '%f %f %f\n', MCS_inv_lon(i), MCS_inv_lat(i), MCS_inv_depth(i));
end

fclose(fid);
