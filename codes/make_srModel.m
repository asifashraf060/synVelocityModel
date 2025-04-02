%% Script for making a srModel based on imported interfaces (basement/moho) 
% load the vp structure (2D--x&z) for upper crust
% load the interface structure (2D--x&y)
% matrix dimension in x direction for interface must be bigger than the vp structure of sediments
% extent of interfaces must be bigger than the vp structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Introduction:
% This script is designed to automate the creation of a seismic velocity model (srModel)
% by integrating multiple geophysical datasets. It leverages imported 2D upper crustal
% velocity matrices and 3D interface structures (basement/moho) 
% to generate a consistent and continuous model grid.
%
% Key features of the script include:
% 1. Data Integration: Loads essential input files such as sediment and upper crust 
%    velocity matrices, interface geometries, station/event locations, elevation data, 
%    and control parameters.
%
% 2. Interpolation & Gap Filling: Employs robust interpolation methods to fill gaps 
%    between the sedimentary data and the interface extents, ensuring that all datasets 
%    align within a unified spatial framework.
%
% 3. Model Assembly: Combines the interpolated data to construct the full 3D velocity 
%    model. This process includes assigning crustal and mantle velocity values, with 
%    optional modifications such as soft-boundary blurring, Vp patching, Gaussian filtering,
%    and checkerboard tests.
%
% 4. Quality Control & Visualization: Conducts thorough NaN checks and produces diagnostic 
%    plots to validate the spatial extent and integrity of the model before final deployment.
%
% 5. Model Saving: Once validated, the completed srModel is saved for use in further seismic 
%    imaging and analysis tasks, particularly within the stingray environment.
%
% Author: Asif (originally created on Apr 20, 2022; updated through 2024-2025)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%clear all, close all, clc

% find input and output dirs
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir)
cd ..
cd outputs/
out_dir      = [pwd, '/entireModel/'];
in_dir_slab  = [pwd, '/slab/structures/'];
in_dir_UC    = [pwd, '/upCrust/structures/'];
cd ..
cd inputs/
in_dir_SR    = [pwd, '/stingray_structures/'];
cd(scriptDir)

addpath(out_dir); addpath(in_dir_slab); addpath(in_dir_UC); addpath(in_dir_SR);

%find the path
theUpCr    = which('upCrust_31-Mar-2025_UpCrust_Vp_12points_highSltzTerrane_AUG13_2024.xlsx.mat');   % 2D vp matrix for upper crust
theInt     = which('int_3D_31-Mar-2025_casieGP_bloch.mat'); % Interfaces
theStation = which('srStation_29-Oct-2024_CG&SG_.mat');
theEvent   = which('srEvent_28-Oct-2024_CG&SG_.mat');
theElevation = which('srElevation_feb2_CGSG.mat');
theGeometry  = which('OR2012_srGeometry.mat');
theControl   = which('srControl_casExp.mat');

load(theUpCr)    %loaded structure name is upCr_2Dmat
load(theInt)     %loaded structure name is int_3Dmat
load(theStation) %loaded structure name is srStation
load(theEvent)   %loaded structure name is srEvent
load(theGeometry)

% Assign paramters to modify standard model-
mantle_vp          = 8.1;
basement_vp_top    = 6.1;
basement_vp_bottom = 6.8;
top_depth          = 0;

xy_spacing         = .5; % input: km
z_spacing          = .5; % input: km
onlyMoho           = 1;  % input: logical
onlyUpCr           = 0;  % input: logical
onlyMantle         = 1;  % input: logical

soft_boundary      = 0;  % input: logical
        % blurring inputs
            window = 10; % windowSize

vp_patch           = 0;  % input: logical
        % Vp Patch input
             vp_ip = 5.5;  

gauss_filter       = 0;  % input: logical
        % gaussian filtering inputs
              g_mg = 20; % magnitude of filtering

checker            = 0;  % input: logical
        % checkerboard test inputs
              x_sp = 40; 
              y_sp = 40; 
              z_sp = 40; % xyz spacing
              mag  = .2; % magnitude of change in percent
              condition = 'whole';

save_model      = 1;

%% AUTO CALCULATION -- no input below this line
%% fill the gaps between interface and sedimentary vp matrix

% fix the model name
model_name = append('srModel_', date ,'_', unq_model_name);

if soft_boundary == 1
    model_name = append(model_name, '_softB_',string(window));
end

if vp_patch == 1
    model_name = append(model_name, '_vPatch',string(vp_ip));
end

if gauss_filter == 1
    model_name = append(model_name, '_gaussFilt', string(g_mg));
end

if checker    == 1
    model_name = append(model_name, '_chck', string(x_sp), string(y_sp), string(x_sp),'M' ,string(mag), condition);
end

disp('Developing Velocity Model ...')

model_depth = min(upCr_2Dmat.zPos);

% create x positions to fill up 2D spaces
xfill_west = min(int_3Dmat.xPos):0.2:min(upCr_2Dmat.xPos);
xfill_east = max(upCr_2Dmat.xPos):0.2:max(int_3Dmat.xPos);
xCombo     = (horzcat(xfill_west, upCr_2Dmat.xPos, xfill_east))';
xfilled    = (min(xCombo):xy_spacing:max(xCombo))';
yfilled    = min(int_3Dmat.yPos):xy_spacing:max(int_3Dmat.yPos);

% fill the spaces for vp

depth4interp     = flip((model_depth):z_spacing: top_depth);
vpfill_west      = repmat(upCr_2Dmat.vp(:,1), size(xfill_west));
vpfill_east      = repmat(upCr_2Dmat.vp(:,end), size(xfill_east));
vp_xCombo        = horzcat(vpfill_west, upCr_2Dmat.vp, vpfill_east);
[xc_grd, zc_grd] = meshgrid(upCr_2Dmat.zPos, xCombo);
[xf_grd, zf_grd] = meshgrid(depth4interp, xfilled);
vp_xfilled       = griddata(xc_grd, zc_grd, vp_xCombo', xf_grd, zf_grd);
vp_xyfilled      = permute(repmat(vp_xfilled, 1, 1, length(yfilled)), [1 3 2]);

disp('1. Gaps been filled to match the imported interface extent')

% make the interface and sedimentary matrix in the same position in 3D space
disp('2. Starting 2-D interpolation ...')

basement_interp     = (interp2(int_3Dmat.xPos, int_3Dmat.yPos', int_3Dmat.basement,...
                               xfilled, yfilled))';
moho_interp         = (interp2(int_3Dmat.xPos, int_3Dmat.yPos', int_3Dmat.moho,...
                               xfilled, yfilled))';
basement_interp_elv = (interp2(int_3Dmat.xPos, int_3Dmat.yPos', int_3Dmat.basement_elv,...
                               xfilled, yfilled))';
moho_interp_elv     = (interp2(int_3Dmat.xPos, int_3Dmat.yPos', int_3Dmat.moho_elv,...
                               xfilled, yfilled))';
elv_interp          = (interp2(int_3Dmat.xPos, int_3Dmat.yPos', int_3Dmat.elevation,...
                               xfilled, yfilled))';

disp('      Finished interpolation')

%% NaN test
disp('3. NaN test ...')
if ~isempty(find(isnan(vp_xyfilled)))
    error('         vp has nan')
end
if ~isempty(find(isnan(basement_interp)))
    error('         basement interface has nan')
end
if ~isempty(find(isnan(moho_interp)))
   error('          moho interface has nan')
end
disp('      There is no nan in the data')

%% fill the crust
disp('4. Filling other tectonic structures ...')
for i = 1:length(yfilled)
    vp_2D_loop    = squeeze(vp_xyfilled (:, i ,:));
    basement_loop = basement_interp(:,i);
    moho_loop     = moho_interp(:,i);
        for j = 1:length(xfilled)
            vp_1D            = vp_2D_loop(j,:);
            basement_loop2   = basement_loop(j);
            moho_loop2       = moho_loop(j);
            cut_index        = find(depth4interp<basement_loop2 & depth4interp>moho_loop2);
            vp_crust         = linspace(basement_vp_top, basement_vp_bottom, length(cut_index));
            vp_1D(cut_index) = vp_crust;
            vp_crust2D(j,:)  = vp_1D;
        end
    vp_crust3D(:,i,:) = vp_crust2D;
end
disp('      Subducting crust been filled')

%% fill the mantle

if onlyMantle == 0 % including or not including the subducting crust in vp model
    vp_input = vp_crust3D;  % vp model with subducting crust
else
    vp_input = vp_xyfilled; % vp model with no subducting crust
end

for m = 1:length(yfilled)
        vp_2D_loop    = squeeze(vp_input (:,m,:));
        basement_loop = basement_interp(:,m);
        moho_loop     = moho_interp(:,m);
             for n = 1:length(xfilled)
                 vp_1D           = vp_2D_loop(n,:);
                 basement_loop2  = basement_loop(n);
                 moho_loop2      = moho_loop(n);
                 cut_index       = find(depth4interp<moho_loop2);
                    if isempty(cut_index)
                        error('deepest moho is deeper than assigned highest moho depth')
                    end
                 vp_moho         = zeros(length(cut_index))+mantle_vp;
                 vp_moho         = vp_moho(1,:);
                 vp_1D(cut_index)= vp_moho;
                 vp_moho2D(n,:)  = vp_1D;
             end
        vp_moho3D(:,m,:) = vp_moho2D;
end
disp('      Upper Mantle been filled')

%% NaN test
disp('5. NaN test ...')
if onlyUpCr == 0
    if ~isempty(find(isnan(vp_moho3D)))
        error('         vp has nan')
    end
end

if onlyUpCr == 1
    if ~isempty(find(isnan(vp_xyfilled)))
        error('         vp has nan')
    end
end
disp('      Vp model has no nan')

%% make the srModel structure
clear srModel

% slowness
if onlyUpCr == 0
    srModel.P.u = 1./vp_moho3D;
else
    srModel.P.u = 1./vp_xyfilled;
end

srModel.xg  = xfilled;
srModel.yg  = yfilled;
srModel.zg  = depth4interp';

%now load srGeometry to convert x and y into lon, lat
[xg_grd, yg_grd] = meshgrid(srModel.xg, srModel.yg);
[srModel.LON, srModel.LAT] = xy2map(xg_grd, yg_grd, srGeometry);

frc_lon = srModel.LON;
frc_lat = srModel.LAT;

if onlyMoho == 1
    %now assign interfaces -- only MOHO
    srModel.interface(1).id                         = 1;
    srModel.interface(1).name                       = {'moho'};
    [srModel.interface(1).Y,srModel.interface(1).X] = meshgrid(srModel.yg, srModel.xg);
    srModel.interface(1).Z                          = moho_interp;
    srModel.interface(1).elevation                  = moho_interp_elv; 
else
%now assign interfaces -- both BASEMENT and MOHO
    srModel.interface(1).id                         = 1;
    srModel.interface(1).name                       = {'basement'};
    [srModel.interface(1).Y,srModel.interface(1).X] = meshgrid(srModel.yg, srModel.xg);
    srModel.interface(1).Z                          = basement_interp;
    srModel.interface(1).elevation                  = basement_interp_elv;
    srModel.interface(2).id                         = 2;
    srModel.interface(2).name                       = {'moho'};
    [srModel.interface(2).Y,srModel.interface(2).X] = meshgrid(srModel.yg, srModel.xg);
    srModel.interface(2).Z                          = moho_interp;
    srModel.interface(2).elevation                  = moho_interp_elv;
end

%now, write ghead
srModel.ghead(1) = srModel.xg(1);
srModel.ghead(2) = srModel.yg(1);
srModel.ghead(3) = length(srModel.xg);
srModel.ghead(4) = length(srModel.yg);
srModel.ghead(5) = length(srModel.zg);
srModel.ghead(6) = abs(srModel.xg(1)-srModel.xg(2));
srModel.ghead(7) = abs(srModel.yg(1)-srModel.yg(2));
srModel.ghead(8) = abs(srModel.zg(1)-srModel.zg(2));

%% PLOTS to check
figure(33), clf
plot(srModel.LON, srModel.LAT)
hold on
plot(srStation.longitude, srStation.latitude)
title('Check Extent--before loading srModel')
saveas(gcf, [out_dir, '/plots/model_extent_b4Load.jpg'])

figure(2), clf
[x_grd, y_grd] = meshgrid(srModel.xg, srModel.zg);
[x, y] = contourf(x_grd, y_grd, (squeeze(1./srModel.P.u (:,(round(length(srModel.yg)/2)),:)))');
colormap(jet)
colorbar
hold on
title('before loading srModel', 'FontSize', 16)
saveas(gcf, [out_dir, '/plots/srModel_b4Load.jpg'])

close all

%% use stingray environment to cross check srModel

disp('Loading through stingray environment...')
% save and load the srModel
theModel = [out_dir, '/models/' ,model_name, '.mat'];
save(theModel, 'srModel')
% load the structures through stingray
srElevation  = load_srElevation(theElevation, srGeometry);
srControl  = load_srControl(theControl);
srModel   = load_srModel(theModel, srControl, srGeometry, srElevation);

% FORCING -- ignoring problems inside load_srModel
srModel.elevation = elv_interp;

% modify the model

if soft_boundary == 1
    disp('6. Blurring the model by convoluting with moving average ...')
    srModel = apply_srModel_mAvg(srModel, window, 'UpDown');
end

if vp_patch == 1

    disp('8. Applying a Vp patch ...')
    disp('      select points on figure')

    srModel.P.u = apply_srModel_vPatch(srModel, vp_ip);
end

if gauss_filter == 1
    disp('7. Applying Gaussian filter ...')

    v = 1./srModel.P.u;
    
    v_filt = imgaussfilt(v, g_mg);
    
    srModel.P.u = 1./v_filt;
end

if checker == 1
    disp('9. Converting into a checkerboard model ...')

    waveDimension = [x_sp y_sp z_sp];
    prcntDv       = mag ;

    srModel.P.u = apply_srModel_checkerboard(srModel, prcntDv, waveDimension, condition, vp_xyfilled);
end

%% NaN test
disp('Final Testing for nans...')
if ~isempty(find(isnan(srModel.elevation)))
    error('srModel.elevation has NaN')
end

if ~isempty(find(isnan(srModel.interface(1).elevation)))
    error('srModel.interface(1).elevation has NaN')
end

if ~isempty(find(isnan(srModel.interface(2).elevation)))
    figure(112391), clf
    
    error('srModel.interface(2).elevation has NaN')
end
disp('      There are no nans')

%% plots to check

figure(22), clf
[x_grd, y_grd] = meshgrid(srModel.xg, srModel.yg);
[x, y] = contourf(x_grd, y_grd, (squeeze(1./srModel.P.u (:,:,(round(length(srModel.zg)/2)))))');
colormap(jet)
colorbar
hold on
title('after loading srModel', 'FontSize', 16)
axis equal
saveas(gcf, [out_dir, '/plots/srModel_afLoad.jpg'])

close all

figure(3), clf
plot(srModel.LON, srModel.LAT, 'y')
hold on
plot(srStation.longitude, srStation.latitude, 'vr')
plot(srEvent.longitude, srEvent.latitude, '.r')
title('Check Extent--after loading srModel')
saveas(gcf, [out_dir, '/plots/model_extent_afLoad.jpg'])

figure(5), clf
contourf(srModel.LON, srModel.LAT, srModel.elevation)
title('check elevation')
saveas(gcf, [out_dir, '/plots/elevation.jpg'])

close all

%% SAVING
if save_model == 1
    save(theModel, 'srModel', '-v7.3')
    disp('MODEL SUCCESSFULLY MADE AND SAVED!')
    disp(theModel)
end