%% Script for developing a smooth 2-D upper crust velocity model
% develops a smooth and denser upper crustal velocity model from few points
% asif, Nov 23, 2022

%% Input

vp_tb  = table2array(readtable('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/UpCrust_Vp_12points_highSltzTerrane_AUG13_2024.xlsx'));
mltpl_factor = 20;
filt_value   = 15;
out_dir      = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/';
model_name   = 'upCrust_2024Aug13_2HighSltz_12points_3kline';
% srGeometry
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/OR2012_srGeometry.mat')

%% Extract info from table

disp('Extracting info from table...')
lon_rw = vp_tb(1,:);    %row for LON values from imported table
lon    = lon_rw(2:end); %ignoring the first value of table
dp_cl  = vp_tb(:,1);    %column for DEPTH values from imported table
dp     = dp_cl(2:end);  %%ignoring the first value of table
vp_cl  = vp_tb(:,(2:end));
vp     = vp_cl((2:end),:);

%% Check longitude and depth values are consistent with imported table
if ~length(vp)==length(dp)
    error('Number of rows in imported table should match the number of depth values')
end
if ~width(vp)==length(lon)
    error('Number of columns in imported table should match the number of longitude values')
end
%% Prepare matrices for scatteredInterpolant function
% all data must be in column vector format
a = [];
for i = 1:width(vp)
    vp_lp = vp(:,i);
    dp_lp = dp;
    ln_lp = repmat(lon(i), length(dp_lp),1);
        if isempty(a)
            a = linspace(1,length(dp_lp), length(dp));
        end
        if ~isempty(a)
            a = a + length(dp_lp);
        end
    vp_scInt(a,1) = vp_lp;
    dp_scInt(a,1) = dp_lp;
    ln_scInt(a,1) = ln_lp;
end
%% Scattered Interpolation
disp('Applying scattered interpolation through Vp matrix...')
f                        = scatteredInterpolant(ln_scInt, dp_scInt, vp_scInt);
lon_finer                = linspace(lon(1), lon(end), length(lon)*mltpl_factor);
dp_matching              = linspace(dp(1), dp(end), length(lon_finer));
[ln_grd_Int, dp_grd_Int] = meshgrid(lon_finer, dp_matching);
vp_Int                   = f(ln_grd_Int, dp_grd_Int);
vp_Int_filt              = imgaussfilt(vp_Int, filt_value);
% plot to check the interpolation
figure(1), clf
[CC,h] = contourf(ln_grd_Int, dp_grd_Int, vp_Int_filt, [min(vp_Int_filt(:)):.1:max(vp_Int_filt(:))]);
colormap(flip(jet))
xlabel('Depth')
ylabel('Longitude')
colorbar
title('Check the interpolated Vp distribution', 'FontSize', 16)
%% Make SED structure
% making a fake lat array to apply map2xy function
fake_lat                 = linspace(40, 45, length(lon_finer));
[lon_fn_grd, lat_fk_grd] = meshgrid(lon_finer, fake_lat);
[xg, yg]                 = map2xy(lon_fn_grd, lat_fk_grd, srGeometry);
xg                       = xg(end,:);
yg                       = yg(:,end);
% developing the structure
upCr_2Dmat.vp   = vp_Int_filt;
upCr_2Dmat.xPos = xg;
upCr_2Dmat.zPos = dp_matching;
upCr_2Dmat.lon  = lon_finer;
save([out_dir, model_name, '.mat'], 'upCr_2Dmat')
disp('upCr_2Dmat structure has been saved')