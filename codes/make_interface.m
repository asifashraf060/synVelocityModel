
close all, clear all, clc

%% Input

% directories of the literature based slab matlab files
% longitude grid
theLon = '/Users/asifashraf/Documents/code_repo/MohoSlabCascadia/outputs/models/LON_grid.mat';
% latitude grid
theLat = '/Users/asifashraf/Documents/code_repo/MohoSlabCascadia/outputs/models/LAT_grid.mat';
% interface grid
theInt = '/Users/asifashraf/Documents/code_repo/MohoSlabCascadia/outputs/models/casiePF_bloch.mat';

% srElevation
theElevation = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srElevation_feb2_CGSG.mat';

% Set the longitude (lnlm) and latitude (ltlm) boundaries for the final merged model.
lon_lim = [-126.5 -122.4]; lat_lim = [41.5 44.5];
lat_lon_interval = .01;


%% Calculation

% find output dir
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir)
cd ..
cd outputs/slab/
out_dir = [pwd, '/'];
cd(scriptDir)

% make arrays for set longitude and latitude values
lon_arr = min(lon_lim(:)):lat_lon_interval:max(lon_lim(:));
lat_arr = min(lat_lim(:)):lat_lon_interval:max(lat_lim(:));
% convert lon-lat into local X-Y
map2xy(lon_arr)
% make grid for lat-lon arrays
[lonG, latG] = meshgrid(lon_arr, lat_arr);

disp('Loading imported slab ...')
% load imported slab
lonS = load(theLon); lonS = lonS.lnAr_g;
latS = load(theLat); latS = latS.ltAr_g;
slab = load(theInt); slab = slab.(string(fieldnames(slab)));

disp('Checking set lat-lon limit ...')
if min(lonS(:))>min(lon_lim) || max(lonS(:))<max(lon_lim)
    error('The set longitude limit falls outside of imported slab!')
end
if min(latS(:))>min(lat_lim) || max(latS(:))<max(lat_lim)
    error('The set latitude limit falls outside of imported slab!')
end

disp('Cutting slab into set lon-lat limit ...')
slabInterp = interp2(lonS, latS, slab, lonG, latG);

% take the elevation out from the slab
disp('Taking elevation out of slab structure')
% load elevation structure
load(theElevation);
elv_lon = [srElevation.header(1), srElevation.header(2)]; elv_lonN = srElevation.header(7);
elv_lat = [srElevation.header(3), srElevation.header(4)]; elv_latN = srElevation.header(8);
elv_lon_arr = linspace(min(elv_lon), max(elv_lon), elv_lonN);
elv_lat_arr = linspace(min(elv_lat), max(elv_lat), elv_latN);
[elv_lnG, elv_ltG] = meshgrid(elv_lon_arr, elv_lat_arr);

% extract elevation in the exported slab grid spacing
elv_extSlab = interp2(elv_lnG, elv_ltG, srElevation.data', lonG, latG);

% subtract the elevation
slab_z = slabInterp - elv_extSlab./1000;

disp('      - Saving plots ...')
% plotting applied elevation
figure(1), clf
contourf(elv_lnG, elv_ltG, srElevation.data')
axis equal
colorbar
saveas(gcf, [out_dir, '/plots/applied_elevation.jpg'])
close all

% plotting imported and exported slab
figure('Renderer', 'painters', 'Position', [20 20 1200 600])

subplot(1,3,1)
contourf(lonS, latS, slab)
hold on
plot([min(lon_lim), max(lon_lim)], [min(lat_lim), min(lat_lim)], '-r', 'LineWidth',2)
plot([min(lon_lim), max(lon_lim)], [max(lat_lim), max(lat_lim)], '-r', 'LineWidth',2)
plot([min(lon_lim), min(lon_lim)], [min(lat_lim), max(lat_lim)], '-r', 'LineWidth',2)
plot([max(lon_lim), max(lon_lim)], [min(lat_lim), max(lat_lim)], '-r', 'LineWidth',2)
colorbar
title('Imported slab')
set(gca, 'FontSize', 16)

subplot(1,3,2)
contourf(lonG, latG, slabInterp)
colorbar
title('Exported slab')
set(gca, 'FontSize', 16)

subplot(1,3,3)
contourf(lonG, latG, slab_z)
colorbar
title('Exported slab (Z)')
set(gca, 'FontSize', 16)

saveas(gcf, [out_dir, '/plots/import_export_slab.jpg'])

close all

% make basement from Moho





