%% CUT interface in x-y plane and add it to another interface

clear all, close all, clc

%% INPUT

% load the int_3Dmat structure, you want to cut and modify
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/output_int_and_UpCr_struc/int_3D_03-Oct-2023int3D_26-Jul-2023int3D_pd16to18_ps01B_td1617_mcsDpCnv_thckTstd_actual.mat');

% STINGRAY INPUT
% srGeomtry
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/OR2012_srGeometry.mat')
% srEvent
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srEvent_20-Feb-2023_PD16to18_PS01Bsub1sub2_D1617.mat');
% srStation
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srStation_21-May-2023_l2to3_g1to6_obsPD16to17rlc.mat')
% path for the slab in any table format (e.g., .txt/.dat/.xls)
theSlab          = '/Users/asifashraf/Dropbox (University of Oregon)/Cascadia_velocity_models/models/mccrory_morph/MCslab_cut_Nov2.txt';
% Elevation data in ascii raster format
theElevation     = '/Users/asifashraf/Documents/Casc_Exp_Files/Elevation/grid123/GMRTv4_1_20221126topo.asc';
% Name the structure
out_dir          = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/output_int_and_UpCr_struc';
unq_outname      = 'mcs_JOIN_mcc_int1';

%% AUTO CAlCULATION

%% DRAW the CUT LINE

disp('Make a cut line on the figure...')

xPos = int_3Dmat.xPos; yPos = int_3Dmat.yPos;

% we need to get the lat lon from x, y position
% if the length of xPos and yPos is not same, xy2map will not work
% that's why to make xPos and yPos length same,
if ~isequal(length(xPos), length(yPos))
    yPos_new = (linspace(yPos(1), yPos(end), length(xPos)))';
    xPos_new = xPos;
end

% lon and lat of the imported interface
[lon_impint, lat_impint]         = xy2map(xPos_new, yPos_new, srGeometry); %xy2map is a function from stingray utilities
[lon_impint_grd, lat_impint_grd] = meshgrid(lon_impint, lat_impint);

% figure to make that cutting line
figure(1), clf
plot(lon_impint_grd, lat_impint_grd)
hold on
plot(srEvent.longitude, srEvent.latitude, '.r')
plot(srStation.longitude, srStation.latitude, 'vr')
set(gca, 'FontSize', 16)
title('Choose several points to make a cutting line')
[lon_cut, lat_cut] = ginput();
plot(lon_cut, lat_cut, '-k', 'LineWidth', 1)
title('Cutting line is made')

%% MAKE lons_cut that cover whole map

disp('Joining the lon value to ends ')

ind_latCut = [];
for i = 1:length(lat_cut)
    ind_latCut(:,i) = find(abs(lat_impint-lat_cut(i)) == min(abs(lat_impint-lat_cut(i))));
end

if ~isequal(sort(ind_latCut), ind_latCut)
    [ind_sorted, I] = sort(ind_latCut);
    ind1            = horzcat(1, ind_sorted , length(lat_impint));
    lon_cut_sorted  = [];
    for i = 1:length(I)
    lon_cut_sorted(i,:) = lon_cut(I(i));
    end
else
    ind1 = horzcat(1, (ind_latCut), length(lat_impint));
end

% make the lon_cut_sorted have same length as ind_sorted
lon_cut_sorted = vertcat(lon_cut_sorted(1), lon_cut_sorted, lon_cut_sorted(end));

lons_cut = [];
Ls = [];
for i = 1:length(ind1)-1
        L = ind1(i):1:ind1(i+1); % length
        M = linspace(lon_cut_sorted(i), lon_cut_sorted(i+1), length(L)-1); % main array
        lons_cut = horzcat(lons_cut, M);
end
lons_cut = horzcat(lons_cut, lons_cut(end));

figure(1)
hold on
plot(lons_cut, lat_impint, '.b', 'LineWidth', 2)

%% extrapolate the basement into new position

[x_posGrd, y_posGrd] = meshgrid(xPos, yPos);
[lon_xPos, lat_yPos] = xy2map(x_posGrd, y_posGrd, srGeometry);

F = scatteredInterpolant(lon_xPos(:), lat_yPos(:), int_3Dmat.basement_elv(:));
basement_new = F(lon_impint_grd, lat_impint_grd);

%% extrapolate the McCrory into new position

slab_tb      = table2array(readtable(theSlab));
mc_lon                    = slab_tb(:,1);
mc_lat                    = slab_tb(:,2);
mc_slab                   = slab_tb(:,3);

F = scatteredInterpolant(mc_lon, mc_lat, mc_slab);
mc_slab_grd = F(lon_impint_grd, lat_impint_grd);

% plot the extrapolated surface
figure(2), clf
plot3(lon_impint_grd, lat_impint_grd, basement_new)
grid on
hold on
plot3(lon_impint_grd, lat_impint_grd, mc_slab_grd, '-y')
plot(lons_cut, lat_impint, '-b', 'LineWidth', 2)

disp('Making the modified interface ...')

basement_modified = [];
for i = 1:length(mc_slab_grd)
    
    bn = basement_new(i,:);
    mc = mc_slab_grd(i,:);
    
    lon_ind = find(abs(lons_cut(i) - lon_impint) == min(abs(lons_cut(i) - lon_impint)));
    
    bn_cut = bn(lon_ind:end); mc_cut = mc(1:(lon_ind-1));
    
    bn_md  = horzcat(mc_cut, bn_cut);
    
    bn_filt = imgaussfilt(bn_md, 5);
    
    basement_modified(i,:) = bn_filt;
    
end

figure(2), hold on
plot3(lon_impint_grd, lat_impint_grd, basement_modified, '-g')


%% Elevation--make Interface.Z

disp('interpolating elevation data...')

[data, metadata] = readgeoraster(theElevation, 'OutputType', 'double');
xmin = metadata.XWorldLimits(1);
xmax = metadata.XWorldLimits(2);
ymin = metadata.YWorldLimits(1);
ymax = metadata.YWorldLimits(2);
xinc = metadata.CellExtentInWorldX;
yinc = metadata.CellExtentInWorldY;
[ny, nx] = size(data);
elv_lon  = [(xmin+xinc):xinc:xmax];
elv_lat  = ([(ymin+yinc):yinc:ymax])';
elv_grd  = flipud(data);
[elv_lon_grd, elv_lat_grd] = meshgrid(elv_lon, elv_lat);
[elv_x_grd, elv_y_grd]     = map2xy(elv_lon_grd, elv_lat_grd, srGeometry);

elevation_interp           = griddata(elv_lon_grd,elv_lat_grd, elv_grd, ... 
                                                lon_impint_grd, lat_impint_grd);
elevation    = elevation_interp/1000; 

figure(5), clf
contourf(lon_impint_grd, lat_impint_grd, (elevation))

%check if basement_Z has any point shallower than zero and make it zero
basement_modified(find(basement_modified > 0)) = -0.01;

%% Make the int_3Dmat structure
clear int_3Dmat

int_3Dmat.basement_elv = basement_modified;
int_3Dmat.basement     = basement_modified - elevation;

int_3Dmat.moho_elv     = int_3Dmat.basement_elv - 6;
int_3Dmat.moho         = int_3Dmat.basement - 6;

int_3Dmat.xPos         = xPos_new;
int_3Dmat.yPos         = yPos_new;

int_3Dmat.elevation    = elevation;

figure(3), clf
plot3(lon_impint_grd, lat_impint_grd, int_3Dmat.basement_elv)
hold on
plot3(lon_impint_grd, lat_impint_grd, int_3Dmat.basement, '-y')
plot3(mc_lon, mc_lat, mc_slab, 'ob')
plot3(lon_impint_grd, lat_impint_grd, int_3Dmat.elevation, '-g')
grid on
plot(srStation.longitude, srStation.latitude, 'vr')
plot(srEvent.longitude, srEvent.latitude, '.r')
ylim([42.5, 43.75])

figure(4), clf
plot3(lon_impint_grd, lat_impint_grd, int_3Dmat.moho_elv)
hold on
plot3(lon_impint_grd, lat_impint_grd, int_3Dmat.moho, '-y')
plot3(mc_lon, mc_lat, mc_slab, 'ob')
plot3(lon_impint_grd, lat_impint_grd, int_3Dmat.elevation, '-g')
grid on
plot(srStation.longitude, srStation.latitude, 'vr')
plot(srEvent.longitude, srEvent.latitude, '.r')
ylim([42.5, 43.75])

figure(5), clf
contourf(lon_impint_grd, lat_impint_grd, int_3Dmat.elevation)
colormap(gray)
colorbar
hold on
plot(srEvent.longitude, srEvent.latitude, '.r')
plot(srStation.longitude, srStation.latitude, 'vy')
plot(lons_cut, lat_impint, '.g', 'LineWidth', 2)
set(gca, 'FontSize', 16)
axis equal
ylim([42.5, 44])




%% Saving int_3Dmat
save(append(out_dir, '/', 'int3D_', date, unq_outname, '.mat'),'int_3Dmat')
disp('~interface structure has been saved~')



























