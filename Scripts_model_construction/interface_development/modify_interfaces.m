% Script to modify interfaces and make four different ones
% modified interfaces will have diffrent depths
% Input is an existing int_3Dmat structure
% 'pct1' and 'pct2' are two percentage inputs that control the depth modification
% 'pct1' is the starting position, left/west of which will be unmodified (typically constrained with some data)
% 'pct2' is the ending position of transition and starting postion of depth
% apply a low filter number (input is 'flt') to smooth corners of transition

clear all, close all, clc

%% INPUT -- MANUAL

% load event and stations
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srStation_13-Jul-2023_l2to3_g1to8_obsPD16to17rlc_pnsn.mat');
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srEvent_20-Feb-2023_PD16to18_PS01Bsub1sub2_D1617.mat');
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/OR2012_srGeometry.mat');

% structure name - int_3Dmat
home_dir = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/output_int_and_UpCr_struc/'; % type a slash at the end
int_str  = 'int3D_16-Nov-2023mcs_bstck_thcktstd';  % type it without .mat

% array knob to set depth change [start_depth end_depth]
inverse_dip = [.5 5];
dip1        = [.5 5];
dip2        = [1 10];
dip3        = [1.5 15];

%number for seperation point from west
pct1 = 30; % percent cut of the model..it will be the starting position
pct2 = 33;
flt  = .5;

%% INPUT -- AUTO
theInt     = append(home_dir, int_str, '.mat');
load(theInt)
%name the surfaces
act_name   = append(int_str, '_actual');
up_name    = append(int_str, '_UP_', string(inverse_dip(1)), 'to', string(inverse_dip(2)));
down_name1 = append(int_str, '_DOWN_', string(dip1(1)), 'to', string(dip1(2)));
down_name2 = append(int_str, '_DOWN_', string(dip2(1)), 'to', string(dip2(2)));
down_name3 = append(int_str, '_DOWN_', string(dip3(1)), 'to', string(dip3(2)));

int = int_3Dmat;
b3D = int.basement + int.elevation;
[X, Y]  = meshgrid(int_3Dmat.yPos, int_3Dmat.xPos);

%Seperation points
sp1  = round(length(int_3Dmat.xPos)*(pct1/100));
sp2  = round(length(int_3Dmat.xPos)*(pct2/100));
mv   = length((sp1+1):1:(sp2-1)); %middle value numbers

%% Loop to create different interfaces
for i = 1:length(int_3Dmat.yPos)
    b1D = flip(b3D(i,:));
    %seperate the slab in X axis
    b_1st  = b1D(1:sp1);
    b_last = b1D(sp2:length(b1D));    
    b_last = b_last';
    
    %make the dip 1 interface
    subt         = ((linspace(dip1(1), dip1(end), length(b_last)))');
    b_last_d     = b_last - subt;
    b_mid_d      = (linspace(b_1st(end), b_last_d(1), mv));
    b_down       = flip(vertcat(b_1st', b_mid_d', b_last_d));
    b_down_filt  = (imgaussfilt(b_down, flt));
    
    %Make more dip 2 interface
    subt_d2      = ((linspace(dip2(1), dip2(end), length(b_last)))');
    b_last_d2    = b_last - subt_d2;
    b_mid_d2     = (linspace(b_1st(end), b_last_d2(1), mv));
    b_down2      = flip(vertcat(b_1st', b_mid_d2', b_last_d2));
    b_down_filt2 = imgaussfilt(b_down2, flt);    

    %Make more dip 3 interface
    subt_d3      = ((linspace(dip3(1), dip3(end), length(b_last)))');
    b_last_d3    = b_last - subt_d3;
    b_mid_d3     = (linspace(b_1st(end), b_last_d3(1), mv));
    b_down3      = flip(vertcat(b_1st', b_mid_d3', b_last_d3(:,1)));
    b_down_filt3 = imgaussfilt(b_down3, flt);  
        
    %make the less dip interface invdip1
    subt_u      = ((linspace(inverse_dip(1), inverse_dip(end), length(b_last)))');
    b_last_u    = b_last + (subt_u);
    b_mid_u     = (linspace(b_1st(end), b_last_u(1), mv));
    b_up        = flip(vertcat(b_1st', b_mid_u', b_last_u(:,1)));
    b_up_filt   = imgaussfilt(b_up, flt);
    
    %store different surfaces
    downline(i,:)  = b_down_filt;
    downline2(i,:) = b_down_filt2;
    downline3(i,:) = b_down_filt3;
    upline(i,:)    = b_up_filt;
    
end

%% figure to check

figure(1), clf
[XX, YY]  = meshgrid(int_3Dmat.xPos, int_3Dmat.yPos);
plot(XX, YY)
hold on
plot(srStation.x, srStation.y, 'vr')
plot(srEvent.x, srEvent.y, '.r')
title('Choose two points to do a profile check')
[x, y] = ginput(2);
plot(x, y, '-k')

xL = linspace(x(1), x(2), 500);
yL = linspace(y(1), y(2), 500);


evL = interp2(XX, YY, int_3Dmat.elevation, xL, yL);


bL  = interp2(XX,YY, b3D, xL, yL);
uL  = interp2(XX,YY, upline, xL, yL);
d1L = interp2(XX,YY, downline, xL, yL);
d2L = interp2(XX,YY, downline2, xL, yL);
d3L = interp2(XX,YY, downline3, xL, yL);

[lonL, latL] = xy2map(xL, yL, srGeometry);



f = figure(22); clf
f.Position = [15 10 500 500]; 
plot(lonL , uL, '-m', 'LineWidth', 2)
hold on
plot(lonL , d1L, '-r', 'LineWidth', 2)
plot(lonL , d2L, '-g', 'LineWidth', 2)
plot(lonL , d3L, '-b', 'LineWidth', 2)
plot(lonL , bL, '-k', 'LineWidth', 2)
plot(lonL, evL, '--k', 'LineWidth', 1)
grid on
set(gca, 'FontSize', 16)

figure(1), clf
plot3(X, Y, downline, 'r')
hold on
grid on
plot3(X, Y, upline, 'b')
plot3(X, Y, downline2, 'g')
plot3(X, Y, downline3, 'c')
plot3(X, Y, b3D)
title('red-dip1 || blue-invdip1 || brown-McCrory || green-dip2 || cyan-dip3')

%% Make the actual interface
%% INT0

int_3Dmat.basement_elv   = b3D;
int_3Dmat.basement       = b3D - int.elevation;
int_3Dmat.moho_elv       = int_3Dmat.basement_elv-6;
int_3Dmat.moho           = int_3Dmat.moho_elv - int.elevation;
name0 = (append(home_dir, 'int_3D_', date, act_name, '.mat'));

save(name0, 'int_3Dmat')

%% MAKE the SHALLOWER surface
%% INVDIP1

int_3Dmat.basement_elv   = upline;
int_3Dmat.basement       = upline - int.elevation;
int_3Dmat.moho_elv       = int_3Dmat.basement_elv-6;
int_3Dmat.moho           = int_3Dmat.moho_elv - int.elevation;
name1 = (append(home_dir, 'int_3D_', date, up_name, '.mat'));

figure(2), clf
plot3(X, Y, int_3Dmat.basement, 'b')
hold on
grid on
plot3(X, Y, int_3Dmat.basement_elv, 'r')
plot3(X, Y, int_3Dmat.moho, 'g')
plot3(X, Y, int_3Dmat.moho_elv, 'y')
title('invdip1, b-base.Z||r-base.ELV||g-moho.Z||y-moho.ELV')

save(name1, 'int_3Dmat')

%% MAKE the DIPPER surface
%% DIP1

int_3Dmat.basement_elv  = downline;
int_3Dmat.basement      = downline - int.elevation;
int_3Dmat.moho_elv      = int_3Dmat.basement_elv-6;
int_3Dmat.moho          = int_3Dmat.moho_elv - int.elevation;
name2 = (append(home_dir, 'int_3D_', date, down_name1, '.mat'));

figure(3), clf
plot3(X, Y, int_3Dmat.basement, 'b')
hold on
grid on
plot3(X, Y, int_3Dmat.basement_elv, 'r')
plot3(X, Y, int_3Dmat.moho, 'g')
plot3(X, Y, int_3Dmat.moho_elv, 'y')
title('dip1, b-base.Z||r-base.ELV||g-moho.Z||y-moho.ELV')

save(name2, 'int_3Dmat')

%% DIP2
int_3Dmat.basement_elv  = downline2;
int_3Dmat.basement      = downline2 - int.elevation;
int_3Dmat.moho_elv      = int_3Dmat.basement_elv-6;
int_3Dmat.moho          = int_3Dmat.moho_elv - int.elevation;
name3 = (append(home_dir, 'int_3D_', date, down_name2, '.mat'));

figure(4), clf
plot3(X, Y, int_3Dmat.basement, 'b')
hold on
grid on
plot3(X, Y, int_3Dmat.basement_elv, 'r')
plot3(X, Y, int_3Dmat.moho, 'g')
plot3(X, Y, int_3Dmat.moho_elv, 'y')
title('dip2, b-base.Z||r-base.ELV||g-moho.Z||y-moho.ELV')

save(name3, 'int_3Dmat')

%% DIP3
int_3Dmat.basement_elv  = downline3;
int_3Dmat.basement      = downline3 - int.elevation;
int_3Dmat.moho_elv      = int_3Dmat.basement_elv-6;
int_3Dmat.moho          = int_3Dmat.moho_elv - int.elevation;
name4 = (append(home_dir, 'int_3D_', date, down_name3, '.mat'));

figure(5), clf
plot3(X, Y, int_3Dmat.basement, 'b')
hold on
grid on
plot3(X, Y, int_3Dmat.basement_elv, 'r')
plot3(X, Y, int_3Dmat.moho, 'g')
plot3(X, Y, int_3Dmat.moho_elv, 'y')
title('dip3, b-base.Z||r-base.ELV||g-moho.Z||y-moho.ELV')

save(name4, 'int_3Dmat')

%% SAVE the interfaces


% 
% 
% figure(3), clf
% plot3(X, Y, int_3Dmat.basement, 'b')
% hold on
% grid on
% plot3(X, Y, int_3Dmat.basement_elv, 'r')
% plot3(X, Y, int_3Dmat.moho, 'g')
% plot3(X, Y, int_3Dmat.moho_elv, 'y')
% title('more dip')
% 
% figure(123), clf
% plot3(srModel.LON, srModel.LAT, srModel.interface(1).elevation, 'b')
% hold on
% plot3(srModel.LON, srModel.LAT, srModel.interface(2).elevation, 'r')


% 
% b1D = b3D(:,1);
% 
% figure(1), clf
% plot(int.xPos, b1D)
% 
% %seperate the slab in X axis
% b_1st  = b1D(1:550);
% b_last = b1D(700:length(b1D));
% 
% %make the dipper interface
% subt   = linspace(2, 30, length(b_last));
% b_last_d = b_last - subt;
% 
% b_mid_d = flip(linspace(b_last_d(1), b_1st(end), (700-551)));
% 
% b_down = vertcat(b_1st, b_mid_d', (b_last_d(1,:))');
% b_down_filt = imgaussfilt(b_down, 50);
% 
% figure(1)
% hold on
% plot(int.xPos, b_down_filt, '-r')
% 
% 
% %make the less dip interface
% subt_u   = linspace(2, 8, length(b_last));
% b_last_u = b_last' + (subt_u);
% b_mid_u  = flip(linspace(b_last_u(1), b_1st(end), (700-551)));
% b_up     = vertcat(b_1st, b_mid_u', b_last_u');
% b_up_filt = imgaussfilt(b_up, 50);
% 
% figure(1)
% hold on
% plot(int.xPos, b_up_filt, '-g')
