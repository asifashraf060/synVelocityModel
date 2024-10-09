%% Script for making interface for imported interface as txt/xlx file
%% Preferred Format of all txt/xlx file - Longitude||Latitude||Depth/TWTT
%% Also can input MCS data to constrain any part
%% This script is using griddata, so it will take time to interpolate
%% Script is prepared for Cascadia subduction zone where lon is decreasing from east to west(min value)
%%                                                   and lat is decreasing from north to south(min value)

clear all, close all, clc

%% INPUT
%% LAT LON for preferred interface extent | Lt = Lat, Ln = Lon
%                                       | n = north, s = south
%                                       | e = east, w = west
Lt_n = 44.48;
Lt_s = 41.81;
Ln_e = -122.29;
Ln_w = -126.52;

saving = 0; % logical to save the interface structure

%% MCS reflection data
theMCS           = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/interface_depthConv_manipulated.txt';
%% STINGRAY INPUT

%srGeomtry
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/OR2012_srGeometry.mat')

%srEvent
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srEvent_28-Apr-2024_SG1234.mat');

%srStation
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srStation16-Apr-2024_SG12&3.mat')

%path for the slab in any table format (e.g., .txt/.dat/.xls) [LON || LAT || DEPTH(-)]

%theSlab          = '/Users/asifashraf/Downloads/data/control-points.txt';

theSlab          = '/Users/asifashraf/Dropbox (University of Oregon)/Cascadia_velocity_models/models/mccrory_morph/MCslab_cut_Apr28_2024.txt';

%Elevation data in ascii raster format
theElevation     = '/Users/asifashraf/Documents/Casc_Exp_Files/Elevation/grid123/GMRTv4_1_20221126topo.asc';
%Name the structure
out_dir          = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/output_int_and_UpCr_struc';
unq_outname      = 'apr28_southExp';

%% CONSISTENT INPUT (change if necessary)
%interface interpolation spacing in km
int_interp_sp = 1;
%column for longitude values in the table
lon_col   = 1;
%column for latitude values in the table
lat_col   = 2;
%column for depth values of the slab
depth_col = 3;

%% AUTO CALCULATION
slab_tb      = (readtable(theSlab));
MCS_oc_tb    = table2array(readtable(theMCS));

%% Calculate the spatial position of the slab

%import the slab structure
mc_lon                    = table2array(slab_tb(:,lon_col));
mc_lat                    = table2array(slab_tb(:,lat_col));
mc_slab                   = (table2array(slab_tb(:,depth_col)));

% filter the spatial extent of slab structure
disp('filtering the imported interface from literature ...')
lt_inp = [Lt_n Lt_s]; ln_inp = [Ln_e Ln_w];
ind_A = find(mc_lon>min(ln_inp)); ind_B = find(mc_lon<max(ln_inp)); ind_ln = intersect(ind_A, ind_B);
ind_A = find(mc_lat>min(lt_inp)); ind_B = find(mc_lat<max(lt_inp)); ind_lt = intersect(ind_A, ind_B);
ind_cut = intersect(ind_ln, ind_lt);
mc_lon = mc_lon(ind_cut); mc_lat = mc_lat(ind_cut); mc_slab = mc_slab(ind_cut);

disp('Calculating x,y for imported interface from literature ...')
[mc_x, mc_y]              = map2xy(mc_lon, mc_lat, srGeometry);
[mc_grd_x, mc_grd_y]      = meshgrid(mc_x, mc_y);
mc_slab_grd               = meshgrid(mc_slab);

%calculate input extent
disp('Calculating x,y based on input lat lons....')
inp_lon                   = [Ln_e, Ln_w];
inp_lat                   = [Lt_n, Lt_s];
[inp_x, inp_y]            = map2xy(inp_lon, inp_lat, srGeometry);
[inp_grd_x, inp_grd_y]    = meshgrid(inp_x, inp_y); 

%x,y for shots and stations
[st_grd_x, st_grd_y]      = map2xy(srStation.longitude, srStation.latitude, srGeometry);
[sh_grd_x, sh_grd_y]      = map2xy(srEvent.longitude, srEvent.latitude, srGeometry);

%PLOT to check if everything is in perferrable range
figure(1), clf
plot(mc_grd_x, mc_grd_y, '.y')
hold on
plot(inp_grd_x, inp_grd_y, '*-g')
title('yellow-imported interface || green-impoted coordinates || red-shots&stations')
hold on
plot(st_grd_x, st_grd_y, 'vr')
plot(sh_grd_x, sh_grd_y, '.r')

%Make grid from input coordinates
inp_x_all                          = flip([inp_x(end):int_interp_sp:inp_x(1)]);
inp_y_all                          = flip([inp_y(end):int_interp_sp:inp_y(1)]);
[inp_x_all_grd, inp_y_all_grd]     = meshgrid(inp_x_all, inp_y_all);
[inp_lon_all_grd, inp_lat_all_grd] = xy2map(inp_x_all, inp_y_all', srGeometry);

%PLOT to check new model excompasses all shots and stations
figure(2), clf
plot(inp_x_all_grd, inp_y_all_grd, 'y')
hold on
plot(st_grd_x, st_grd_y, 'vr')
plot(sh_grd_x, sh_grd_y, '.r')
title('yellow = new model extent')


%% Interpolate into all x and y positions
F                              = scatteredInterpolant(mc_x, mc_y, mc_slab);
[inp_x_all_grd, inp_y_all_grd] = meshgrid(inp_x_all, inp_y_all);
F.Method                       = 'linear';
MCslab_interp                  = F(inp_x_all_grd, inp_y_all_grd);

% forcing filtering
MCslab_interp = imgaussfilt(MCslab_interp, 20);

%plot the interpolated basement
figure(4), clf
plot3(inp_x_all_grd, inp_y_all_grd, MCslab_interp, '.')
hold on
plot3(mc_x, mc_y, mc_slab, 'ob')
grid on
hold on
plot(st_grd_x, st_grd_y, 'vr')
plot(sh_grd_x, sh_grd_y, '.r')
title('New interpolated slab with ScatteredInterpolant', 'FontSize', 16)


%% APPLYING MCS
disp('applying the imported MCS data...')
% imported MCS data
MCS_lon   = MCS_oc_tb(:,1);
MCS_lat   = MCS_oc_tb(:,2);
MCS_depth = MCS_oc_tb(:,3);


% add with bathymetry
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
elv_interp1 = griddata(elv_lon_grd, elv_lat_grd, elv_grd, MCS_lon, MCS_lat); elv_interp1 = elv_interp1/1000;

figure(45), clf
plot3(MCS_lon, MCS_lat, MCS_depth, '.')
grid on
[MCS_x, MCS_y]   = map2xy(MCS_lon, MCS_lat, srGeometry);

% defining extent of interpolation
MCS_lon_interp   = min(MCS_lon-.1):.02:max(MCS_lon+.01);
MCS_lat_interp   = min(MCS_lat-.1):.02:max(MCS_lat+.1);

%% Adding a fake LINES to fill data space of the imported MCS

disp('Adding fake lines ...')

% fake NORTH line
disp('   ... in North')
ind1              = find(MCS_lat == max(MCS_lat));
fakeLineN_lat     = MCS_lat(min(ind1));
fakeLineN_lat_lon = MCS_lon(min(ind1));
ind2              = find(MCS_lon == min(MCS_lon));
fakeLineN_lon     = MCS_lon(min(ind2));
fakeLineN_depth_W = MCS_depth(min(ind2));
fakeLineN_depth_E = MCS_depth(min(ind1));

fakeLineN_lons    = (fakeLineN_lon:.01:fakeLineN_lat_lon)';
fakeLineN_lats    = zeros(length(fakeLineN_lons),1)+ fakeLineN_lat;
fakeLineN_depths  = (linspace(fakeLineN_depth_W, fakeLineN_depth_E, length(fakeLineN_lons)))';

% fake SOUTH line (fsl)
disp('   ... in South')
fsl_lon = min(MCS_lon):.01:max(MCS_lon);
fsl_lat = zeros(length(fsl_lon)) + min(MCS_lat_interp); fsl_lat = fsl_lat(1,:);

% loop to go through each lat lon and assign a depth value
fsl_dp = [];
for i = 1:(length(fsl_lon)-1)
    ln1  = fsl_lon(i); ln2 = fsl_lon(i+1);
    ind1 = find(MCS_lon>=ln1 & MCS_lon<=ln2);
    
    lt3  = MCS_lat(min(ind1));
    dp1  = MCS_depth(min(ind1));
    
    ind2 = find(lt3<max(MCS_lat) & lt3>min(MCS_lat));
    dp2  = dp1(min(ind2));
    
    fsl_dp(:,i)  =  mean(dp2);
    if i == length(fsl_lon)-1
        fsl_dp(:,i+1) = mean(dp2);
    end
end
fsl_dp = fillmissing(fsl_dp, 'nearest'); % to replace the nan values

% fake east line (fel)
disp('   ... in East')
fel_lat = (min(MCS_lat):.01:max(MCS_lat));
fel_lon = zeros(length(fel_lat))+max(MCS_lon_interp); fel_lon = fel_lon(1,:);

% loop to go through each lat lon and assign a depth value
fel_dp = [];
for i = 1:(length(fel_lat)-1)
    lt1  = fel_lat(i); lt2 = fel_lat(i+1);
    ind1 = find(MCS_lat>=lt1 & MCS_lat<=lt2);
    
    lt3  = MCS_lon(min(ind1));
    dp1  = MCS_depth(min(ind1));
    
    ind2 = find(lt3<max(MCS_lon) & lt3>max(MCS_lon)-0.1);
    dp2  = dp1(min(ind2));
    
    fel_dp(:,i) = mean(dp2);
    if i == length(fel_lat)-1
        fel_dp(:,i+1) = mean(dp2);
    end
end
fel_dp = fillmissing(fel_dp, 'nearest'); % to replace the nan values


% fake WEST line (fwl)
disp('   ... in West')
fwl_lat = (min(MCS_lat):.01:max(MCS_lat));
fwl_lon = zeros(length(fel_lat))+min(MCS_lon_interp); fwl_lon = fwl_lon(1,:);

% loop to go through each lat lon and assign a depth value
fwl_dp = [];
for i = 1:(length(fwl_lat)-1)
    lt1  = fwl_lat(i); lt2 = fwl_lat(i+1);
    ind1 = find(MCS_lat>=lt1 & MCS_lat<=lt2);
    
    lt3  = MCS_lon(min(ind1));
    dp1  = MCS_depth(min(ind1));
    
    ind2 = find(lt3>min(MCS_lon) & lt3<min(MCS_lon)+.1);
    dp2  = dp1(min(ind2));
    
    fwl_dp(:,i) = mean(dp2);
    if i == length(fwl_lat)-1
        fwl_dp(:,i+1) = mean(dp2);
    end
end
fwl_dp = fillmissing(fwl_dp, 'nearest'); % to replace the nan values


%% Combine

% adding all the data togther
MCS_lon_ext      = vertcat(MCS_lon, fakeLineN_lons, fsl_lon', fel_lon', fwl_lon');
MCS_lat_ext      = vertcat(MCS_lat, fakeLineN_lats, fsl_lat', fel_lat', fwl_lat');
MCS_depth_ext    = vertcat(MCS_depth, fakeLineN_depths, fsl_dp', fel_dp', fwl_dp');

disp('interpolating over imported MCS data...')

[MCS_lon_interp_grd, MCS_lat_interp_grd] = meshgrid(MCS_lon_interp, MCS_lat_interp);

clear F
F                = scatteredInterpolant(MCS_lon_ext, MCS_lat_ext, MCS_depth_ext);
F.Method         = 'linear';
MCS_interp       = F(MCS_lon_interp_grd, MCS_lat_interp_grd);

MCS_interp_filt  = imgaussfilt(MCS_interp,.5);


% plot
figure(47), clf

            plot3(MCS_lon_ext, MCS_lat_ext, MCS_depth_ext, '.r')
            grid on
            hold on
            plot3(MCS_lon_interp_grd, MCS_lat_interp_grd, MCS_interp, '.b')
            title('Unfiltered interface')

figure(447), clf
            plot3(MCS_lon_ext, MCS_lat_ext, MCS_depth_ext, '.r')
            hold on
            grid on
            plot3(MCS_lon_interp_grd, MCS_lat_interp_grd, MCS_interp_filt, '.g')
            title('filtered interface')

%% Fill MCS dataset to match the input coordinates
%  lon is decreasing from east to west(min value)
%  lat is decreasing from north to south(min value)

disp('Filling MCS dataset to match the input coordinates...')

%  step 1: fill the WEST side
Fln_wMCS                     = (Ln_w-.01):.02:min(unique(MCS_lon_interp_grd)); % FIlling the longitudes in west of MCS grid
Flt_wMCS                     = unique(MCS_lat_interp_grd);               % FIlling the latitudes in west of MCS grid
Fint_wMCS                    = MCS_interp_filt(:,(find((unique(MCS_lon_interp_grd)...
                                                          == min(unique(MCS_lon_interp_grd)))))); % Finding the interface values to fill
Fint_wMCS_grd                = repmat(Fint_wMCS, [1 length(Fln_wMCS)]);
Fint_wMCS_total              = horzcat(Fint_wMCS_grd, MCS_interp_filt);  %combining west fill and MCS grid
[Fln_wMCS_grd, Flt_wMCS_grd] = meshgrid(Fln_wMCS, Flt_wMCS);

figure(665), clf

            plot3(MCS_lon_interp_grd, MCS_lat_interp_grd, MCS_interp_filt, '.g')

            hold on

            plot3(Fln_wMCS_grd, Flt_wMCS_grd, Fint_wMCS_grd, 'b')

            title('Check the filled part from MCS in blue and red', 'FontSize', 16)
            grid on

            xlabel('longitude')
            ylabel('latitude')

%  step 2: fill the NORTH side
Fln_nMCS                     = horzcat(Fln_wMCS, (unique(MCS_lon_interp_grd))'); % FIlling the longitudes in north of MCS grid
Flt_nMCS                     = (max(unique(MCS_lat_interp_grd)):0.02:(Lt_n+.1))';     % FIlling the latitudes in north of MCS grid
Fint_nMCS                    = Fint_wMCS_total(find((unique(MCS_lat_interp_grd) == max(unique(MCS_lat_interp_grd)))),:);
Fint_nMCS_grd                = repmat(Fint_nMCS, [length(Flt_nMCS) 1]);
[Fln_nMCS_grd, Flt_nMCS_grd] = meshgrid(Fln_nMCS, Flt_nMCS);
figure(665)
hold on
plot3(Fln_nMCS_grd, Flt_nMCS_grd, Fint_nMCS_grd, 'r')
%  step 3: fill the SOUTH side
Fln_sMCS                     = horzcat(Fln_wMCS, (unique(MCS_lon_interp_grd))'); % FIlling the longitudes in south of MCS grid
Flt_sMCS                     = ((Lt_s-.1):0.02:min(unique(MCS_lat_interp_grd)))';     % FIlling the latitudes in south of MCS grid
Fint_sMCS                    = Fint_wMCS_total(find((unique(MCS_lat_interp_grd) == min(unique(MCS_lat_interp_grd)))),:);
Fint_sMCS_grd                = repmat(Fint_sMCS, [length(Flt_sMCS) 1]);
[Fln_sMCS_grd, Flt_sMCS_grd] = meshgrid(Fln_sMCS, Flt_sMCS);
figure(665)
hold on
plot3(Fln_sMCS_grd, Flt_sMCS_grd, Fint_sMCS_grd, 'r')
% step 4: COMBINING all the filled parts
MCS_f_ln                     = Fln_nMCS;
MCS_f_lt                     = vertcat(Flt_nMCS, Flt_wMCS, Flt_sMCS);
MCS_f_int                    = vertcat(Fint_nMCS_grd, Fint_wMCS_total, Fint_sMCS_grd);
[MCS_f_ln_grd, MCS_f_lt_grd] = meshgrid(MCS_f_ln, MCS_f_lt);
% step 5: INTERPOLATE into input lat-lon space
interp_lon                       = MCS_f_ln;
interp_lat                       = inp_lat_all_grd(:,1);
[interp_lon_grd, interp_lat_grd] = meshgrid(interp_lon,interp_lat);
MCS_f_int_interp                 = griddata(MCS_f_ln_grd, MCS_f_lt_grd, MCS_f_int, ...
                                                    interp_lon_grd, interp_lat_grd);

figure(666), clf
            plot3(interp_lon_grd, interp_lat_grd, MCS_f_int_interp, 'b.')
            title('Check the MCS data with filled parts combined', 'FontSize', 16)
            grid on
%% Cut the MC interface in accordance with the MCS data
disp('Cutting the imported interface...')
ind            = find(inp_lon_all_grd(end,:)>(max(max(MCS_f_ln_grd))+.5));
MCslab_cut     = MCslab_interp(:,ind);
MCslab_cut_lon = inp_lon_all_grd(:,ind);
MCslab_cut_lat = inp_lat_all_grd(:,ind);

%% Join the MCS and cut MC slab
disp('Joining MCS and McCrory data...')
%easternmost MCS depth
MCS_eastEnd_lon = max(interp_lon);
MCS_eastEnd     = MCS_f_int_interp(:,(find(interp_lon == MCS_eastEnd_lon)));
%westernmost MC depth
MCcut_lon_all   = (MCslab_cut_lon(1,:));
MC_westEnd_lon  = min(MCcut_lon_all);
MC_westEnd      = MCslab_cut(:,(find(MCcut_lon_all == MC_westEnd_lon)));
join_lon        = MCS_eastEnd_lon:.02:MC_westEnd_lon;
clear join_int
for i = 1:length(MC_westEnd)
    MCS = MCS_eastEnd(i);
    MC  = MC_westEnd(i);
    int = linspace(MCS, MC, length(join_lon));
    join_int(i,:) = int;
end
[join_lon_grd, join_lat_grd] = meshgrid(join_lon, interp_lat);

%% Combining all the developed interfaces
disp('Combining all developed interfaces...')
combo_lon       = horzcat(interp_lon_grd, join_lon_grd, MCslab_cut_lon);
combo_lat       = horzcat(interp_lat_grd, join_lat_grd, MCslab_cut_lat);
combo_int       = horzcat(MCS_f_int_interp, join_int, MCslab_cut);

% scatter interpolation
F               = scatteredInterpolant(combo_lon(:), combo_lat(:), combo_int(:));
F.Method        = 'linear';
final_int       = F(inp_lon_all_grd, inp_lat_all_grd);
final_int_filt  = imgaussfilt(final_int, 1);
diff_int_filt   = final_int-final_int_filt;


[combo_x, combo_y] = map2xy(inp_lon_all_grd, inp_lat_all_grd, srGeometry);

figure(669), clf
            plot3(inp_lon_all_grd, inp_lat_all_grd, final_int_filt, '.')
            grid on
            title('filtered interface || green - McCrory || blue - seismic reflection || red - stations & events', 'FontSize', 16)
            hold on
            plot3(mc_lon, mc_lat, mc_slab, 'og')
            grid on
            hold on
            plot(srStation.longitude, srStation.latitude, '*r')
            plot(srEvent.longitude, srEvent.latitude, '.r')
            plot3(MCS_lon, MCS_lat, MCS_depth, '.b')
            xlabel('longitude')
            ylabel('latitude')
            zlabel('elevation')

figure(770), clf
            [CC, h] = contourf(inp_lon_all_grd, inp_lat_all_grd,diff_int_filt);
            clabel(CC)
            colormap('jet')
            colorbar
            title('diff betn interface and filtered interface', 'FontSize', 16)
            hold on
            plot(srStation.longitude, srStation.latitude, '*r')
            hold on
            plot(srEvent.longitude, srEvent.latitude, '.r')

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
elevation_interp           = griddata(elv_x_grd,elv_y_grd, elv_grd, ... 
                                                inp_x_all_grd, inp_y_all_grd);
figure(1112), clf
            subplot(1,2,1)
            contourf(elv_x_grd, elv_y_grd, elv_grd, [-5000:200:3000])
            colormap(winter)
            cc = colorbar;
            hold on
            plot(st_grd_x, st_grd_y, 'vr')
            plot(sh_grd_x, sh_grd_y, '.r')
            title('imported elevation')

            subplot(1,2,2)
            contourf(inp_x_all_grd, inp_y_all_grd, elevation_interp, [-3500:200:3000])
            colormap(winter)
            cc = colorbar;
            hold on
            plot(st_grd_x, st_grd_y, 'vr')
            plot(sh_grd_x, sh_grd_y, '.r')
            title('imported interpolated elevation')
            
elevation    = (elevation_interp/1000); 
basement_elv = final_int_filt;

% make sure elevation and basement are consistent
b_2D_corrected = [];
for i = 1:length(inp_y_all)
    
    ev_1D = elevation(i,:);
    b_1D  = basement_elv(i,:);
    
    b_1D_corrected = [];
    for j = 1:length(b_1D)
        
        ev_p = ev_1D(j);
        b_p  = b_1D(j);
        
        if b_p>ev_p
            b_p = ev_p;
        end
        
        b_1D_corrected(j) = b_p;
    end
    
    b_2D_corrected(i,:) = b_1D_corrected;
end

basement_elv = b_2D_corrected;

basement_Z   = basement_elv - elevation;

figure(5), clf
            plot3(inp_x_all_grd, inp_y_all_grd, basement_elv, '.')
            hold on
            plot3(inp_x_all_grd, inp_y_all_grd, basement_Z, 'y')
            grid on
            plot3(inp_x_all_grd, inp_y_all_grd, elevation, '.g')
            hold on
            plot(st_grd_x, st_grd_y, 'vr')
            plot(sh_grd_x, sh_grd_y, '.r')
            plot(MCS_x, MCS_y, 'ob')
            plot3(mc_x, mc_y, mc_slab, 'ob')

%check if basement_Z has any point shallower than zero and make it zero
basement_Z(find(basement_Z > 0)) = -0.01;

%% Make int_3Dmat structure
int_3Dmat.xPos          = (inp_x_all');
int_3Dmat.yPos          = (inp_y_all');
int_3Dmat.basement      = basement_Z;
int_3Dmat.basement_elv  = basement_elv;

int_3Dmat.moho          = imgaussfilt((basement_Z), 8) - 6;
int_3Dmat.moho_elv      = imgaussfilt((basement_elv), 8) - 6;
int_3Dmat.elevation     = elevation;


figure(990), clf
            plot(int_3Dmat.xPos, int_3Dmat.basement, '--k')
            title('check if the basement in the right direction')

figure(991), clf
            [int_X, int_Y] = meshgrid(int_3Dmat.xPos, int_3Dmat.yPos);
            [int_LON, int_LAT] = xy2map(int_X, int_Y, srGeometry);
            plot(int_LON, int_LAT, '.')
            hold on
            plot(srStation.longitude, srStation.latitude, 'vr')
            plot(srEvent.longitude, srEvent.latitude, '.r')

figure(992), clf
            plot3(int_LON, int_LAT, int_3Dmat.basement_elv, 'b')
            hold on
            grid on
            plot3(int_LON, int_LAT, int_3Dmat.moho_elv, 'g')
            plot3(MCS_lon, MCS_lat, MCS_depth, '.r')
            plot3(mc_lon, mc_lat, mc_slab, 'ok')
            plot(srStation.longitude, srStation.latitude, 'vm')
            plot(srEvent.longitude, srEvent.latitude, '.m')
            plot3(int_LON, int_LAT, int_3Dmat.elevation, 'y')
            title('Basement [black-McCrory] [red-MCS]')
            set(gca, 'FontSize', 16)
            

figure(993), clf
            plot3(int_LON, int_LAT, int_3Dmat.basement_elv, 'b')
            hold on
            grid on
            plot(srStation.longitude, srStation.latitude, 'vm')
            plot(srEvent.longitude, srEvent.latitude, '.m')
            plot3(int_LON, int_LAT, int_3Dmat.elevation, 'y')
            plot3(MCS_lon_ext, MCS_lat_ext, MCS_depth_ext, '.r')
            title('Basement [blue-elevation] [green-Z]')
            set(gca, 'FontSize', 16)




%% Saving int_3Dmat
if saving == 1
    save(append(out_dir, '/', 'int3D_', date, unq_outname, '.mat'),'int_3Dmat')
    disp('~interface structure has been saved~')
end


