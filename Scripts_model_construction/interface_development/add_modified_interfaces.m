%CONNECT TWO INTERFACES
%imported interfaces must be in the same space
%asif, Jan 22, 2023
clear all, close all
addpath '/Users/asifashraf/Documents/MATLAB/'
addpath '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/';
addpath '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/'

%% INPUT
Int_file1 = 'int_3D_20-Dec-2022int3D_25-Nov-2022MCS_McCrory_l3k2k_g123_DOWN_1to10.mat';  %file name for interface structure in North
Int_file2 = 'int_3D_22-Jan-2023_int0_deeper_mod.mat'; %file name for interface structure in South
out_name  = append('int_3D_', date, '_int0and1_modadd.mat');

%% INPUT -- consistent
Station   = 'srStation_20-Nov-2022s3000+grid1+grid2+grid3+s2000.mat';
Event     = 'srEvent_11-Nov-2022PD16+PS01B+PD18.mat';
Geometry  = 'OR2012_srGeometry.mat';

%% AUTO
int1 = load(which(Int_file1)); %Northern interface
int2 = load(which(Int_file2)); %Southern interface
load(which(Station));
load(which(Event));
load(which(Geometry));

%convert x,y to lat,lon
[int_xy, int_yx] = meshgrid(int1.int_3Dmat.xPos, int1.int_3Dmat.yPos);
[int_lon_grd, int_lat_grd] = xy2map(int_xy, int_yx, srGeometry);
int_lon = (int_lon_grd(round(length(int1.int_3Dmat.xPos)/2),:)'); %taking the middle row
int_lat = int_lat_grd(:,round(length(int1.int_3Dmat.yPos)/2));    %taking the middle column

figure(1),clf
plot(int_lon_grd, int_lat_grd)
hold on
plot(srStation.longitude, srStation.latitude, 'vr')
plot(srEvent.longitude, srEvent.latitude, '.r')
[lon_input, lat_input] = ginput(1);

%find the index for the nearest value to 'lat_input'
d   = abs(int_lat - lat_input);
ind = find(d == min(d));
%find all the indices for cutting int_lat
ind_st_N     = find(int_lat == max(int_lat)); %maximum because the experiment is in North Hemisphere
ind_N_series = [ind_st_N:1:ind];
ind_st_S     = find(int_lat == min(int_lat)); %minimum because the experiment is in North Hemisphere
ind_S_series = [ind+1:1:ind_st_S];

%% Basement
int_north_bZ   = int1.int_3Dmat.basement(ind_N_series,:);
int_south_bZ   = int2.int_3Dmat.basement(ind_S_series,:);
int_add_bZ     = vertcat(int_north_bZ, int_south_bZ);
int_north_bELV = int1.int_3Dmat.basement_elv(ind_N_series,:);
int_south_bELV = int2.int_3Dmat.basement_elv(ind_S_series,:);
int_add_bELV   = vertcat(int_north_bELV, int_south_bELV);
%% Moho
int_north_mZ   = int1.int_3Dmat.moho(ind_N_series,:);
int_south_mZ   = int2.int_3Dmat.moho(ind_S_series,:);
int_add_mZ     = vertcat(int_north_mZ, int_south_mZ);
int_north_mELV = int1.int_3Dmat.moho_elv(ind_N_series,:);
int_south_mELV = int2.int_3Dmat.moho_elv(ind_S_series,:);
int_add_mELV   = vertcat(int_north_mELV, int_south_mELV);
%% add smooth TRANSITION
lt_trn_st  = lat_input+.1;              lt_trn_end  = lat_input-.1;
d1         = abs(int_lat - lt_trn_st);  d2          = abs(int_lat - lt_trn_end);
ind_trn_st = find(d1==min(d1));         ind_trn_end = find(d2==min(d2));
ind_trn_series = ind_trn_st:1:ind_trn_end;
%in basement
clear basement_Z
for i = 1:length(int1.int_3Dmat.xPos)
        b      = int_add_bZ(:,i);
        zV_trn = b(ind_trn_series);
        zV_new = linspace(zV_trn(1), zV_trn(end), length(ind_trn_series));
        b(ind_trn_series) = zV_new;
        basement_Z(:,i) = b;
end
clear basement_ELV
for i = 1:length(int1.int_3Dmat.xPos)
        b      = int_add_bELV(:,i);
        zV_trn = b(ind_trn_series);
        zV_new = linspace(zV_trn(1), zV_trn(end), length(ind_trn_series));
        b(ind_trn_series) = zV_new;
        basement_ELV(:,i) = b;
end
%in moho
clear moho_Z
for i = 1:length(int1.int_3Dmat.xPos)
        b      = int_add_mZ(:,i);
        zV_trn = b(ind_trn_series);
        zV_new = linspace(zV_trn(1), zV_trn(end), length(ind_trn_series));
        b(ind_trn_series) = zV_new;
        moho_Z(:,i) = b;
end
clear moho_ELV
for i = 1:length(int1.int_3Dmat.xPos)
        b      = int_add_mELV(:,i);
        zV_trn = b(ind_trn_series);
        zV_new = linspace(zV_trn(1), zV_trn(end), length(ind_trn_series));
        b(ind_trn_series) = zV_new;
        moho_ELV(:,i) = b;
end
figure(2), clf
plot3(int_lon_grd, int_lat_grd, basement_ELV)
hold on
plot3(int_lon_grd, int_lat_grd, moho_ELV, '-r')
grid on
hold on
plot(srStation.longitude, srStation.latitude, 'vr')
plot(srEvent.longitude, srEvent.latitude, '.r')

%% MAKE NEW INTERFACE STRUCTURE
int_3Dmat.xPos         = int1.int_3Dmat.xPos;
int_3Dmat.yPos         = int1.int_3Dmat.yPos;
int_3Dmat.basement     = basement_Z;
int_3Dmat.basement_elv = basement_ELV;
int_3Dmat.moho         = moho_Z;
int_3Dmat.moho_elv     = moho_ELV;
int_3Dmat.elevation    = int1.int_3Dmat.elevation;
save(append('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/', out_name),'int_3Dmat')

