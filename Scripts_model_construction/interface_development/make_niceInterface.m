%%
clear all, close all, clc

out_dir   = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/';


theMCS    = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/interface_depthConv_manipulated.txt';
MCS_oc_tb = table2array(readtable(theMCS));

% imported MCS data
MCS_lon   = MCS_oc_tb(:,1);
MCS_lat   = MCS_oc_tb(:,2);
MCS_depth = MCS_oc_tb(:,3);


%%

f = scatteredInterpolant(MCS_lon, MCS_lat, MCS_depth);

lat = 42.4:.1:44.2;
lon = -126.2:.1:-124.4;

[lon_grd, lat_grd] = meshgrid(lon, lat);

int_interp = f(lon_grd, lat_grd);

figure(1), clf
plot3(lon_grd(:), lat_grd(:), int_interp(:), '.b')
hold on
grid on
plot3(MCS_lon, MCS_lat, MCS_depth, 'or')

lon = lon_grd(:); lat = lat_grd(:); dp = int_interp(:);



%% Write

fid = fopen(append(append(out_dir, 'interface_nice.txt')), 'w');

if fid == -1
    error('Unable to opne file')
end

for i = 1:length(lon)
    fprintf(fid, '%f %f %f\n', lon(i), lat(i), dp(i));
end

fclose(fid);