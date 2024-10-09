%% MODIFY MOHO INTERFACE FROM srMODEL STRUCTURE
clear, close all, clc
% 231025_112721 - int1
% 231025_112750 - int2
% 231025_112950 - int3
% 231025_113249 - int4
% 231025_113422 - int5 (McCrory)

% 231120_165725 = TOMO1 (BStck)
% 231120_184126 = TOMO2
% 231120_223547 = TOMO3

%% INPUT

% directory of the interface INVERTED srMODEL
inv_model_dir = '/Users/asifashraf/Talapas/tlOutput/231120_165725/srModel_it5.mat';
thePert       = '/Users/asifashraf/Talapas/tlOutput/231120_165725/tlPert_it5.mat';

% directory of the INITIAL srMODEL that was inverted to make the above inverted model
int_model_dir = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/output_new_models/srModel_20-Nov-2023_int_mcs_bstck_addedTogether.mat';

% directory of stingray structures
station       = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srStation_3-Jul-2023_l2to3_g1to6_obsPD16to17rlc_pnsn.mat';
event         = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srEvent_20-Feb-2023_PD16to18_PS01Bsub1sub2_D1617.mat';

% name of the INTERFACE we want to modify (moho/basement/seafloor)
name_int      = 'moho';

int_name      = 'interface 5';

%% CALCULATION -- AUTO

% load the srModel structures as initial (INI) and inverted (INV) interface
inv = load(inv_model_dir);
int = load(int_model_dir);
load(station); load(event);

z = inv.srModel.zg;

% get the index for the moho interface
ind = find(contains(cat(length(inv.srModel.interface), ...
                    inv.srModel.interface.name), name_int) == 1);
                
% extract both inverted and intial interface 
inv_int = inv.srModel.interface(ind); int_int = int.srModel.interface(ind);

% extract lat lon for these interfaces
LON = int.srModel.LON; LAT = int.srModel.LAT;

% make a masked interface
int_diff = (inv_int.elevation - int_int.elevation);
    % loop to make put NaNs on values less than 0.5
    int_diff_masked = [];
    for i = 1:length(int_diff)
         row = (int_diff(i,:));
         ind = find(abs(row)<0.5);
         row(ind) = nan;
        int_diff_masked(i,:) = row; 
    end


% choose points to make a graph for inverted interface
figure(2), clf
[C, h] = contourf(LON, LAT, int_diff_masked);
clabel(C, h)
colorbar
hold on
plot(srStation.longitude, srStation.latitude, 'vr')
plot(srEvent.longitude, srEvent.latitude, '.r')
title('Choose points to normal pressure')
set(gca, 'FontSize', 16)
[lon1, lat1] = ginput;

figure(3), clf
[C, h] = contourf(LON, LAT, int_diff_masked);
clabel(C, h)
colorbar
hold on
plot(srStation.longitude, srStation.latitude, 'vr')
plot(srEvent.longitude, srEvent.latitude, '.r')
title('Choose points to to make an inverted interface')
set(gca, 'FontSize', 16)
[lon2, lat2] = ginput;


% convert that lat lon to x, y
[x1, y1] = map2xy(lon1, lat1, inv.srModel.srGeometry);
[x2, y2] = map2xy(lon2, lat2, inv.srModel.srGeometry);



% using the choosen points extract moho depth from inverted interface
int_inv_points1 = griddata(LON, LAT, inv_int.elevation, lon1, lat1);
int_inv_points2 = griddata(LON, LAT, inv_int.elevation, lon2, lat2);

% basement from moho depth
basement_fromInv = int_inv_points2 + 6;

%% McCrory
% using the choosen points extract moho depth from McCrory interface
load('/Users/asifashraf/Dropbox (University of Oregon)/Cascadia_velocity_models/models/mccrory_morph/MC_Slab.mat');
    % Make the McCrory surface NaN free
    nan_ind = unique(vertcat(find(isnan(MC_Lon) == 1), ... 
                                find(isnan(MC_Lat) == 1), ...
                                    find(isnan(MC_depth) == 1)));
    MC_Lon(nan_ind) = []; MC_Lat(nan_ind) = []; MC_depth(nan_ind) = [];

        % use Scattered Interpolation to extact Mc points on choosen points
        F = scatteredInterpolant(MC_Lon, MC_Lat, MC_depth);
        int_mc_points1  = F(lon2, lat2);

        
figure(339), clf
plot3(lon2, lat2, int_mc_points1, 'ob')
hold on
plot3(lon2, lat2, basement_fromInv, 'or')
grid on
        
lon_spacings = linspace(min(lon2),max(lon2), 50);
lat_spacings = linspace(min(lat2),max(lat2), 50);

[lon_sGrd, lat_sGrd] = meshgrid(lon_spacings, lat_spacings);

% scatteredInterpolation
F = scatteredInterpolant(lon2, lat2, basement_fromInv);
inv_int_sc = F(lon_sGrd, lat_sGrd);
inv_bs_sc_filt = imgaussfilt(inv_int_sc, 2);

F = scatteredInterpolant(lon2, lat2, int_mc_points1);
mc_int_sc  = F(lon_sGrd, lat_sGrd);
mc_int_sc_filt = imgaussfilt(mc_int_sc, 2);



figure(123), clf

plot3(lon_sGrd, lat_sGrd, inv_bs_sc_filt)
grid on
hold on
plot3(lon_sGrd, lat_sGrd, mc_int_sc_filt, '-y')
plot3(lon2, lat2, basement_fromInv, 'or')

hold on

s1  = shaperead('/Users/asifashraf/Documents/GIS_shapeFiles/all_states/States_shapefile/States_shapefile.shp');
mapshow(s1, 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)

s2  = shaperead('/Users/asifashraf/Documents/Casc_Exp_Files/siletz terrane shapefile/In WGS 84/SltzTr_W84.shp');
mapshow(s2, 'EdgeColor', 'c', 'FaceColor', 'g', 'FaceAlpha', 1, 'LineWidth', 2)

xlim([-125 -124]); ylim([42.5 43.5]); zlim([-25 0])


for i = 1:length(lon2)
    v11(i,:) = horzcat(lon2(i), lat2(i), int_inv_points2(i));
end

disp(v11)



% 
% figure(3), clf
% 
% subplot(1, 3, 1)
% plot(lat1, int_inv_points1, '-vr', 'LineWidth', 2)
% hold on
% plot(lat1, (int_mc_points1 - 6), '--vr', 'LineWidth', 1)
% ylim([-25 -10])
% grid on
% title('Inverted (-) || McCrory (--) ')
% set(gca, 'XDir', 'reverse', 'FontSize', 16)
% 
% subplot(1, 3, 2)
% [C, h] = contourf(LON, LAT, int_diff_masked);
% clabel(C, h)
% colorbar
% hold on
% plot(srStation.longitude, srStation.latitude, 'vy')
% plot(srEvent.longitude, srEvent.latitude, '.y')
% title('Difference (inverted and initial)')
% set(gca, 'FontSize', 16)
% plot(lon1, lat1, '-or', 'LineWidth', 4)

%% NORMAL PRESSURE

% to make load DWS
load(thePert)
display('interpolating tlPert dws in model space...')
[PertX, PertY, PertZ]    = meshgrid(tlPert.U.P.y, tlPert.U.P.x, tlPert.U.P.z);
[ModelX, ModelY, ModelZ] = meshgrid(inv.srModel.yg, inv.srModel.xg, inv.srModel.zg);
[ModelXZ, ModelZX]       = meshgrid(inv.srModel.xg, inv.srModel.zg);
[ModelXY, ModelYX]       = meshgrid(inv.srModel.xg, inv.srModel.yg);
dws_modelSpace           = interp3(PertX, PertY, PertZ, tlPert.U.P.dws, ModelX, ModelY, ModelZ);

% Normal pressure calculation
np_3D = [];
for i = 1:length(inv.srModel.yg)
    up_yLine  = squeeze(inv.srModel.P.u(:,i,:));
    int_yLine = squeeze(int.srModel.interface(1).Z(:,i));
    dws_yLine = squeeze(dws_modelSpace(:,i,:));
 
    np_2D = [];
    for j = 1:length(inv.srModel.xg)
        
        dws_xLine = dws_yLine(j,:);
        
        % extact 2-D Vp line
        up_zLine = up_yLine(j,:);
        vp_zLine = 1./up_zLine;
        % extract the basement depth and find the index for it
        int_xPoint = int_yLine(j)+6;
        dp_ind = (find(abs(inv.srModel.zg- int_xPoint) == ...
                                min(abs(inv.srModel.zg- int_xPoint))));
        % use interface depth to choose Vp values of upper crust
        vp_upCrust = vp_zLine(1:dp_ind);
        % convert the up_crust Vp to density and to normal pressure
        np_point = [];
        for k = 1:length(vp_upCrust)
            vp_cnv = vp_upCrust(k); % Vp value subject to conversion
            if vp_cnv < 6
                rho = 2.4372 + 0.0761 * vp_cnv;
            else
                rho = 1.74 * vp_cnv^(0.25);
            end
            np_point(:,k) = rho; % np = normal pressure
        end
        
        dws_sum_array  = cumsum(dws_xLine);
        dws_point      = dws_sum_array(end);
        np_sum_array   = cumsum(np_point);
        
        if dws_point>0
            np_2D(j,:) = np_sum_array(end)*abs(z(1)-z(2))*9.81;
        else
            np_2D(j,:) = nan;
        end
    end
    np_3D(:,i) = np_2D;
end

figure(4), clf
[cc, h] = contourf(inv.srModel.LON, inv.srModel.LAT, np_3D);
clabel(cc)
hold on
plot(lon1, lat1, '*r')


[np_sP] = griddata(inv.srModel.LON, inv.srModel.LAT, np_3D, lon1, lat1);
% 
% 
% figure(3)
% hold on
% subplot(1, 3, 3)
% plot(lat1, np_sP, '-vm', 'LineWidth', 2)
% grid on
% title('Normal Pressure')
% set(gca, 'XDir', 'reverse', 'FontSize', 16)
% 
% figure(5), clf
% 
% subplot(1, 2, 1)
% plot(lat1, int_inv_points1, '-vr', 'LineWidth', 2)
% hold on
% plot(lat1, (int_mc_points1 - 6), '--vr', 'LineWidth', 1)
% ylim([-25 -10])
% grid on
% 
% set(gca, 'XDir', 'reverse', 'FontSize', 16)
% 
% 
% subplot(1, 2, 2)
% plot(lat1, np_sP, '-vm', 'LineWidth', 2)
% grid on
% 
% set(gca, 'XDir', 'reverse', 'FontSize', 16)


%% Mask the inverted and initial interface

dint = abs(inv.srModel.interface(1).elevation - int.srModel.interface(1).elevation);

[row, column, data] = find(dint<2);

inverted_interface   = -1 .* (inv.srModel.interface(1).elevation);
inverted_interface2  = -1 .* (inv.srModel.interface(1).elevation);

initial_interface    = -1 .* (int.srModel.interface(1).elevation);

for i = 1:length(row)
    
    inverted_interface(row(i), column(i)) = nan;
    initial_interface(row(i), column(i)) = nan;

    disp(i)
end


%%
clear h

figure(11), clf
% [c,h] = contour(int.srModel.LON, int.srModel.LAT, ...
%                     inverted_interface2, 'Fill', 'on');
%                 
% hFills = h.FacePrims;
% [hFills.ColorType] = deal('truecoloralpha');            
% AlphaGradient=exp(linspace(log(10),log(100),length(inverted_interface2))); %The actual values that I want to use for the alpha values
% for idx = 1:numel(hFills)
%     hFills(idx).ColorData(4) = AlphaGradient(idx);
% end
% 
% eventFcn = @(srcObj, e) updateTransparency(srcObj);
% addlistener(h, 'MarkedClean', eventFcn);

hold on

[c,h] = contour(int.srModel.LON, int.srModel.LAT, ...
                 inverted_interface2, ...
                    [10:2:30],...
                        '--');

hold on

[cc,h] = contourf(int.srModel.LON, int.srModel.LAT, inverted_interface, [10:2:30]);
%clabel(cc, 'FontSize', 15)
clabel(c, 'FontSize', 15, 'Color', 'r')
clabel(cc, 'FontSize', 15, 'Color', 'r')


colormap(cool); colorbar
hold on
plot(srEvent.longitude, srEvent.latitude, '.g')
plot(srStation.longitude, srStation.latitude, 'vg', 'MarkerFaceColor','g')
xlim([-125.2 -124]); ylim([42.5 43.6])
title('Inverted interface')

hold on

s1  = shaperead('/Users/asifashraf/Documents/GIS_shapeFiles/all_states/States_shapefile/States_shapefile.shp');
mapshow(s1, 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)

s2  = shaperead('/Users/asifashraf/Documents/Casc_Exp_Files/siletz terrane shapefile/In WGS 84/SltzTr_W84.shp');
mapshow(s2, 'EdgeColor', 'c', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)

set(gca, 'FontSize', 16)

%%
figure(111), clf

[c,h] = contour(int.srModel.LON, int.srModel.LAT, ...
                 inverted_interface2, ...
                    [10:2:30],...
                        '--');

hold on

[cc,h] = contourf(int.srModel.LON, int.srModel.LAT, inverted_interface, [10:2:30]);
%clabel(cc, 'FontSize', 15)
%clabel(c, 'FontSize', 15, 'Color', 'r')
%clabel(cc, 'FontSize', 15, 'Color', 'r')


colormap(cool); colorbar
hold on
plot(srEvent.longitude, srEvent.latitude, '.g')
plot(srStation.longitude, srStation.latitude, 'vg', 'MarkerFaceColor','g')
xlim([-125.2 -124]); ylim([42.5 43.6])
title('Inverted interface')

hold on

s1  = shaperead('/Users/asifashraf/Documents/GIS_shapeFiles/all_states/States_shapefile/States_shapefile.shp');
mapshow(s1, 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)

s2  = shaperead('/Users/asifashraf/Documents/Casc_Exp_Files/siletz terrane shapefile/In WGS 84/SltzTr_W84.shp');
mapshow(s2, 'EdgeColor', 'r', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)

set(gca, 'FontSize', 16)


%%
figure(12), clf
[cc,h] = contourf(int.srModel.LON, int.srModel.LAT, initial_interface);
clabel(cc, 'FontSize', 15)
colormap(cool); colorbar
hold on
plot(srEvent.longitude, srEvent.latitude, '.r')
plot(srStation.longitude, srStation.latitude, 'vr')
xlim([-125.2 -124]); ylim([42.5 43.6])
title('Initial interface')

hold on

s1  = shaperead('/Users/asifashraf/Documents/GIS_shapeFiles/all_states/States_shapefile/States_shapefile.shp');
mapshow(s1, 'EdgeColor', 'g', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)

s2  = shaperead('/Users/asifashraf/Documents/Casc_Exp_Files/siletz terrane shapefile/In WGS 84/SltzTr_W84.shp');
mapshow(s2, 'EdgeColor', 'm', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)



set(gca, 'FontSize', 16)

%%
figure(13), clf
diff_int = initial_interface-inverted_interface;
[cc,h]   = contourf(int.srModel.LON, int.srModel.LAT, diff_int, [round(min(min(diff_int))):1:round(max(max(diff_int)))]);
clabel(cc, 'FontSize', 15)
colormap(cool); colorbar
hold on
plot(srEvent.longitude, srEvent.latitude, '.g')
plot(srStation.longitude, srStation.latitude, 'vg', 'MarkerFaceColor','g')
xlim([-125.2 -124]); ylim([42.5 43.6])
title('Initial-Inverted')

hold on

s1  = shaperead('/Users/asifashraf/Documents/GIS_shapeFiles/all_states/States_shapefile/States_shapefile.shp');
mapshow(s1, 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)

s2  = shaperead('/Users/asifashraf/Documents/Casc_Exp_Files/siletz terrane shapefile/In WGS 84/SltzTr_W84.shp');
mapshow(s2, 'EdgeColor', 'c', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)



set(gca, 'FontSize', 16)


disp(append('Results are for', int_name))

%%
figure(133), clf
diff_int = initial_interface-inverted_interface;
[cc,h]   = contourf(int.srModel.LON, int.srModel.LAT, diff_int, [round(min(min(diff_int))):1:round(max(max(diff_int)))]);
%clabel(cc, 'FontSize', 15)
colormap(cool); colorbar
hold on
plot(srEvent.longitude, srEvent.latitude, '.g')
plot(srStation.longitude, srStation.latitude, 'vg', 'MarkerFaceColor','g')
xlim([-125.2 -124]); ylim([42.5 43.6])
title('Initial-Inverted')

hold on

s1  = shaperead('/Users/asifashraf/Documents/GIS_shapeFiles/all_states/States_shapefile/States_shapefile.shp');
mapshow(s1, 'EdgeColor', 'k', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)

s2  = shaperead('/Users/asifashraf/Documents/Casc_Exp_Files/siletz terrane shapefile/In WGS 84/SltzTr_W84.shp');
mapshow(s2, 'EdgeColor', 'r', 'FaceColor', 'w', 'FaceAlpha', .01, 'LineWidth', 2)



set(gca, 'FontSize', 16)


disp(append('Results are for', int_name))


%% DWS Calculation

for i = 1:length(inv.srModel.yg)
    
    dws_2D = squeeze(dws_modelSpace(:,i,:));
    dws_2D(find(dws_2D<10000)) = nan;
    
    
    [xg_grid, zg_grid] = meshgrid(inv.srModel.xg, inv.srModel.zg);
    if ~isempty(find(dws_2D>0))
        figure(110), clf
        contourf(xg_grid, zg_grid, dws_2D')
    else
        
    end
end














