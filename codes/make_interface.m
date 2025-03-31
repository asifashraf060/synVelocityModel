%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   This script processes a geophysical slab model by reading in input 
%   data (including Moho depth, elevation, and crustal geometry) and 
%   generating a detailed interface structure, that can be used to make the
%   3-D synthetic velocity model.
%
% Usage:
%   Update the input file names and set the geographic boundaries by 
%   modifying the variables in the "Input" section. The script requires 
%   the following files (located in the 'inputs/slab/' directory):
%       - casieGP_bloch.mat      : Interface structure file.
%       - LON_grid.mat           : Longitude grid.
%       - LAT_grid.mat           : Latitude grid.
%       - srElevation_feb2_CGSG.mat : Elevation data.
%       - OR2012_srGeometry.mat  : Geometry and mapping parameters.
%
% Output:
%   The script produces:
%       - A merged 3D interface structure (saved in the outputs directory).
%       - Several plots visualizing the processed slab and basement surfaces.
%
% Author: Asif Ashraf
% Date: March 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all, clear all, clc

%% Input

% filename of the interface structure
theInt = 'casieGP_bloch.mat';

% Set the longitude (lnlm) and latitude (ltlm) boundaries for the final merged model.
lon_lim = [-126.5 -122.4]; lat_lim = [41.5 44.5];
lat_lon_interval = .01;

% output structure name
unq_name = 'casieGP_bloch';

%% Calculation

% find output dir
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir)
cd ..
cd outputs/slab/
out_dir = [pwd, '/'];
cd ../..
cd inputs/
in_sl_dir  = [pwd, '/slab/'];
in_sr_dir  = [pwd, '/stingray_structures/'];
cd(scriptDir)

theLon = [in_sl_dir, 'LON_grid.mat']; theLat = [in_sl_dir, 'LAT_grid.mat'];
theInt = [in_sl_dir, theInt];
theElevation = [in_sr_dir, 'srElevation_feb2_CGSG.mat'];
theGeometry  = [in_sr_dir, 'OR2012_srGeometry.mat'];

% make arrays for set longitude and latitude values
lon_arr = min(lon_lim(:)):lat_lon_interval:max(lon_lim(:));
lat_arr = min(lat_lim(:)):lat_lon_interval:max(lat_lim(:));
[lonG, latG]   = meshgrid(lon_arr, lat_arr);
% convert lon-lat into local X-Y
load(theGeometry);
[xG, yG] = map2xy(lonG, latG, srGeometry);

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

% make basement from Moho
[Xb, Yb, Zb]    = surfNormal(xG, yG, slabInterp, 6);
[Lonb, Latb]    = xy2map(Xb, Yb, srGeometry);
[basement]      = griddata(Lonb, Latb, Zb, lonG, latG);

[Xbz, Ybz, Zbz] = surfNormal(xG, yG, slab_z, 6);
[Lonbz, Latbz]  = xy2map(Xbz, Ybz, srGeometry);
[basementZ]     = griddata(Lonbz, Latbz, Zbz, lonG, latG);

basement  = fillmissing(basement, 'linear');
basementZ = fillmissing(basementZ, 'linear');

% combine all grids into int_3Dmat structure

% fake lon and lat array to input into map2xy function
fake_lon = linspace(min(lon_arr), max(lon_arr), length(lat_arr));
fake_lat = linspace(min(lat_arr), max(lat_arr), length(lon_arr));

[xPos, ~] =  map2xy(lon_arr, fake_lat, srGeometry);
[~, yPos] =  map2xy(fake_lon, lat_arr, srGeometry);

int_3Dmat.xPos = xPos';
int_3Dmat.yPos = yPos';
int_3Dmat.basement = basementZ;
int_3Dmat.basement_elv = basement;
int_3Dmat.moho = slab_z;
int_3Dmat.moho_elv = slabInterp;

int_3Dmat.elevation = elv_extSlab./1000;

disp('      - Saving plots ...')
% plotting applied elevation
figure(1), clf
contourf(elv_lnG, elv_ltG, srElevation.data')
axis equal
colorbar
saveas(gcf, [out_dir, '/plots/applied_elevation.jpg'])
close all

% plotting imported and exported slab
figure('Renderer', 'painters', 'Position', [20 20 1400 650])

subplot(1,5,1)
contourf(lonS, latS, slab)
hold on
plot([min(lon_lim), max(lon_lim)], [min(lat_lim), min(lat_lim)], '-r', 'LineWidth',2)
plot([min(lon_lim), max(lon_lim)], [max(lat_lim), max(lat_lim)], '-r', 'LineWidth',2)
plot([min(lon_lim), min(lon_lim)], [min(lat_lim), max(lat_lim)], '-r', 'LineWidth',2)
plot([max(lon_lim), max(lon_lim)], [min(lat_lim), max(lat_lim)], '-r', 'LineWidth',2)
colorbar
title('Imported moho')
set(gca, 'FontSize', 16)

subplot(1,5,2)
[C,h] = contourf(lonG, latG, slabInterp);
clabel(C,h)
colorbar
title('Exp. moho')
set(gca, 'FontSize', 16)

subplot(1,5,3)
[C, h] = contourf(lonG, latG, slab_z);
clabel(C, h)
colorbar
title('Exp. moho (Z)')
set(gca, 'FontSize', 16)

subplot(1,5,4)
[C, h] = contourf(lonG, latG, basement);
clabel(C, h)
colorbar
title('Exp. basement')
set(gca, 'FontSize', 16)

subplot(1,5,5)
[C, h] = contourf(lonG, latG, basementZ);
clabel(C, h)
colorbar
title('Exp. basement (Z)')
set(gca, 'FontSize', 16)

saveas(gcf, [out_dir, '/plots/import_export_slab.jpg'])

close all


figure('Renderer', 'painters', 'Position', [20 20 1400 650])

[xPG, yPG] = meshgrid(int_3Dmat.xPos, int_3Dmat.yPos);

subplot(1,4,1)
[C,h] = contourf(xPG, yPG, int_3Dmat.moho_elv);
clabel(C,h)
colorbar
title('mohoElv')
set(gca, 'FontSize', 16)

subplot(1,4,2)
[C,h] = contourf(xPG, yPG, int_3Dmat.moho);
clabel(C,h)
colorbar
title('moho')
set(gca, 'FontSize', 16)

subplot(1,4,3)
[C,h] = contourf(xPG, yPG, int_3Dmat.basement_elv);
clabel(C,h)
colorbar
title('basementElv')
set(gca, 'FontSize', 16)

subplot(1,4,4)
[C,h] = contourf(xPG, yPG, int_3Dmat.basement);
clabel(C,h)
colorbar
title('basement')
set(gca, 'FontSize', 16)

saveas(gcf, [out_dir, '/plots/exported_structure.jpg'])

close all

disp('      - Saving structure ...')

save([out_dir, '/structures/int_3D_', date, '_' unq_name, '.mat'], 'int_3Dmat')

disp(int_3Dmat)