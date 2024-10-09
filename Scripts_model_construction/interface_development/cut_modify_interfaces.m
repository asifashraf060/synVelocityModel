%MODIFY ONE INTERFACE
%asif, Jan 22 2023
clear all, close all
addpath '/Users/asifashraf/Documents/MATLAB/'
addpath '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/';
addpath '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/'

%% INPUT
Int_file      = 'int_3D_22-Jan-2023_int0_deeper_mod.mat';  %file name for interface structure
ln_cut        = -124.65;   % longitudinal value to cut the interface
cut_direction = 2;         % 1 for low-high, 2 for high-low in the xPos
depth_modify  = 0;         % depth value in km for modifying the interface || Positive Value = modify into deeper interface & vice versa
finer_cut     = 1;         %logical for finer cut
save          = 0;         % logical for save option, 1 = yes; 0 = no
out_name      = append('int_3D_', date, '_int0_deeper_mod.mat');

%% INPUT -- consistent
Station  = 'srStation_20-Nov-2022s3000+grid1+grid2+grid3+s2000.mat';
Event    = 'srEvent_11-Nov-2022PD16+PS01B+PD18.mat';
Geometry = 'OR2012_srGeometry.mat';


%% INPUT -- finer cut
lon_jump1   = -124.693; %must be bigger than lon_jump2
lon_trans1  = -124.644;
lon_jump2   = -124.847;
lon_trans2  = -125.25;
d_array_L2H = [-.5 -1]; %array containing 2 depth values in km that will be applied || Positive Value = modify into deeper interface & vice versa
                      %onto the interface between lon_jump1 & lon_jump2
                      %1st value is for lower longitude value
                      %2nd value is for higher longitude value

%% AUTO
load(which(Int_file));
load(which(Station));
load(which(Event));
load(which(Geometry));
disp('Files are loaded')
%convert x,y to lat,lon
[int_xy, int_yx] = meshgrid(int_3Dmat.xPos, int_3Dmat.yPos);
[int_lon_grd, int_lat_grd] = xy2map(int_xy, int_yx, srGeometry);
int_lon = (int_lon_grd(round(length(int_3Dmat.xPos)/2),:)'); %taking the middle row
int_lat = int_lat_grd(:,round(length(int_3Dmat.yPos)/2));    %taking the middle column
%find the index for the nearest value to 'ln_cut'
d   = abs(int_lon - ln_cut);
ind = find(d == min(d));
%find all the indices for cutting int_lon
if cut_direction == 1
    ind_st = find(int_lon == min(int_lon)); 
else
    ind_st = find(int_lon == max(int_lon));
end
inds   = ind_st:1:ind;
if isempty(inds)
    inds = ind:1:ind_st;
end

disp('Indices are collected')
disp('Cutting the interface...')
%CUT all the ind_3Dmat fields according to the indices
basement_Z_cut     = int_3Dmat.basement(:,inds);
basement_elv_cut   = int_3Dmat.basement_elv(:,inds);
moho_Z_cut         = int_3Dmat.moho(:,inds);
moho_elv_cut       = int_3Dmat.moho_elv(:,inds);
%MODIFY rest of the interface
basement_Z_mdfy    = int_3Dmat.basement(:,setdiff(1:end, inds)) - depth_modify;
basement_elv_mdfy  = int_3Dmat.basement_elv(:,setdiff(1:end, inds)) - depth_modify;
moho_Z_mdfy        = int_3Dmat.moho(:,setdiff(1:end, inds)) - depth_modify;
moho_elv_mdfy      = int_3Dmat.moho_elv(:,setdiff(1:end, inds)) - depth_modify;
%ADD the cut and modified interface
if cut_direction == 1
    basement_Z_new     = horzcat(basement_Z_mdfy, basement_Z_cut);
    basement_elv_new   = horzcat(basement_elv_mdfy, basement_elv_cut);
    moho_Z_new         = horzcat(moho_Z_mdfy, moho_Z_cut);
    moho_elv_new       = horzcat(moho_elv_mdfy, moho_elv_cut);
else
    basement_Z_new     = horzcat(basement_Z_cut, basement_Z_mdfy);
    basement_elv_new   = horzcat(basement_elv_cut, basement_elv_mdfy);
    moho_Z_new         = horzcat(moho_Z_cut, moho_Z_mdfy);
    moho_elv_new       = horzcat(moho_elv_cut, moho_elv_mdfy);  
end
%TRANSITION smoother
ln_cut_trans       = ln_cut - .1; % longitudinal value for to start the transition
d                  = abs(int_lon - ln_cut_trans); %find the index for the nearest value to 'ln_cut_transition'
ind_trans_st       = find(d == min(d));
ln_cut_trans       = ln_cut + .1; % longitudinal value for to start the transition
d                  = abs(int_lon - ln_cut_trans); %find the index for the nearest value to 'ln_cut_transition'
ind_trans_end      = find(d == min(d));
ind_trans_series   = ind_trans_end:1:ind_trans_st;
%% APPLY TRANSITION
clear basement_Z_new_trans
for i = 1:length(int_3Dmat.yPos)
        bZ                          = basement_Z_new(i,:);
        zV_trans                    = bZ(ind_trans_series);
        zV_new                      = linspace(zV_trans(1), zV_trans(end),...
                                                length(ind_trans_series));
        bZ(ind_trans_series)        = zV_new;
        basement_Z_new_trans(i,:)   = bZ;
end
clear basement_elv_new_trans
for i = 1:length(int_3Dmat.yPos)
        bZ                          = basement_elv_new(i,:);
        zV_trans                    = bZ(ind_trans_series);
        zV_new                      = linspace(zV_trans(1), zV_trans(end),...
                                                length(ind_trans_series));
        bZ(ind_trans_series)        = zV_new;
        basement_elv_new_trans(i,:) = bZ;
end
clear moho_Z_new_trans
for i = 1:length(int_3Dmat.yPos)
        bZ                          = moho_Z_new(i,:);
        zV_trans                    = bZ(ind_trans_series);
        zV_new                      = linspace(zV_trans(1), zV_trans(end),...
                                                length(ind_trans_series));
        bZ(ind_trans_series)        = zV_new;
        moho_Z_new_trans(i,:)       = bZ;
end
clear moho_elv_new_trans
for i = 1:length(int_3Dmat.yPos)
        bZ                          = moho_elv_new(i,:);
        zV_trans                    = bZ(ind_trans_series);
        zV_new                      = linspace(zV_trans(1), zV_trans(end),...
                                                length(ind_trans_series));
        bZ(ind_trans_series)        = zV_new;
        moho_elv_new_trans(i,:)     = bZ;
end

figure(1), clf
hold on
plot3(int_lon_grd, int_lat_grd, int_3Dmat.moho_elv, 'g')
title('Green=Previous Interface || Yellow = Modified interface')
grid on
set(gca, 'FontSize', 16)
%% SURGICAL CUT
yPos = int_3Dmat.yPos;
if finer_cut == 1
[basement_Z_new_trans_jumped,...
       basement_elv_new_trans_jumped,...
          moho_Z_new_trans_jumped,...
            moho_elv_new_trans_jumped]... 
                = interface_surgical_cut(int_lon, lon_jump1, lon_jump2, ...
                                            lon_trans1, lon_trans2, d_array_L2H,...
                                                yPos, basement_Z_new_trans,...
                                                    basement_elv_new_trans, moho_Z_new_trans,...
                                                        moho_elv_new_trans);
end

%% MAKE NEW INTERFACE STRUCTURE
if finer_cut == 1
int_3Dmat.basement     = basement_Z_new_trans_jumped;
int_3Dmat.basement_elv = basement_elv_new_trans_jumped;
int_3Dmat.moho         = moho_Z_new_trans_jumped;
int_3Dmat.moho_elv     = moho_elv_new_trans_jumped;
else
int_3Dmat.basement     = basement_Z_new_trans;
int_3Dmat.basement_elv = basement_elv_new_trans;
int_3Dmat.moho         = moho_Z_new_trans;
int_3Dmat.moho_elv     = moho_elv_new_trans;
end

if save == 1
    save(append('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/', out_name),'int_3Dmat');
else
    figure(1)
    hold on
    plot3(int_lon_grd, int_lat_grd, int_3Dmat.moho_elv, 'y')
end


