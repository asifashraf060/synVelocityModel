%% Script for making basement interface from given Moho depth
% This script creates basement interface by applying a constant thickness
% For dipping interface thickness is calculated by applying some geometry

%% INPUT
close all, clear all, clc

% int_3D.mat matrix
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/output_int_and_UpCr_struc/int3D_16-Nov-2023mcs_bstck.mat');

% adjust this array -- KNOB
arr = [2.15 2.15];

%% Calculation

x = int_3Dmat.xPos;
y = int_3Dmat.yPos;

m = int_3Dmat.moho_elv;
b = int_3Dmat.basement_elv;

ind_mid = find(y == y(round(length(y)/2)));

b_line  = b(ind_mid,:); m_line = m(ind_mid,:);

figure(1),clf
plot(x, b_line, '.b')
hold on
plot(x, m_line, '.m')
title('original interfaces')

% check for potential threshold to start the complex thickness calculation
dip = [];
for i = 1:(length(m_line)-1)
    diff = m_line(i) - m_line(i+1);
    dip(:,i) = diff;
end
figure(2), clf
plot(x(1:(length(x)-1)), dip, 'o-b')
gca.YDir = 'reverse';
title('Choose 1 point')
[px, py] = ginput(1);

% find the indices for cut and uncut
ind_xCut = find(abs(x-px-20) == min(abs(x-px-20)));
ind_xUncut = find(abs(x-px+20) == min(abs(x-px+20)));

% modify the cut array
m_mdf               = m_line(1:ind_xCut);
arr_mdf             = -1 * linspace(arr(1), arr(2), length(m_mdf));
m_mdf               = m_mdf+arr_mdf;
m_line(1:ind_xCut)  = m_mdf;

% modify the array between cut and uncut
m_mdf               = m_line(ind_xCut:ind_xUncut);
arr_mdf             = linspace(m_mdf(1), m_mdf(end), length(m_mdf));
m_line(ind_xCut:ind_xUncut)  = arr_mdf;

figure(3),clf
plot(x, b_line, '.b')
hold on
plot(x, m_line, '.m')
title('Modified interfaces')


%% Apply the change to the 3-D interface
new_moho_elv = [];
for i = 1:length(int_3Dmat.yPos)
    
    % extract a 2-D line
    m_line = int_3Dmat.moho_elv(i,:); m_line_z = int_3Dmat.moho(i,:);
    
    % find the indices for cut and uncut
    ind_xCut   = find(abs(x-px-20) == min(abs(x-px-20)));
    ind_xUncut = find(abs(x-px+20) == min(abs(x-px+20)));
    
    % modify the cut array (m_line)
    m_mdf               = m_line(1:ind_xCut);
    arr_mdf             = -1 * linspace(arr(1), arr(2), length(m_mdf));
    m_mdf               = m_mdf+arr_mdf;
    m_line(1:ind_xCut)  = m_mdf;

    % modify the cut array (m_line_z)
    m_mdf_z             = m_line_z(1:ind_xCut);
    arr_mdf_z           = -1 * linspace(arr(1), arr(2), length(m_mdf_z));
    m_mdf_z             = m_mdf_z+arr_mdf_z;
    m_line_z(1:ind_xCut) = m_mdf_z;
    
    % modify the array between cut and uncut (m_line)
    m_mdf               = m_line(ind_xCut:ind_xUncut);
    arr_mdf             = linspace(m_mdf(1), m_mdf(end), length(m_mdf));
    m_line(ind_xCut:ind_xUncut)  = arr_mdf;
    
    % modify the array between cut and uncut (m_line_z)
    m_mdf_z             = m_line_z(ind_xCut:ind_xUncut);
    arr_mdf_z           = linspace(m_mdf_z(1), m_mdf_z(end), length(m_mdf_z));
    m_line_z(ind_xCut:ind_xUncut) = arr_mdf_z;
    
    new_moho_elv(i,:) = m_line;
    new_moho_z(i,:)   = m_line_z;
end

figure(4), clf
[y_grd, x_grd] = meshgrid(int_3Dmat.yPos, int_3Dmat.xPos);
plot3(y_grd, x_grd, new_moho_z)
hold on
plot3(y_grd, x_grd, int_3Dmat.basement, 'y')

%% Save

int_3Dmat.moho     = new_moho_z;
int_3Dmat.moho_elv = new_moho_elv;

save('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/output_int_and_UpCr_struc/int3D_16-Nov-2023mcs_bstck_thcktstd.mat', 'int_3Dmat')

