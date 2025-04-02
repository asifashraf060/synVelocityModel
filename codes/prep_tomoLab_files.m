% Script to prepare all required files to run Stingray & TomoLab inversion


clear all, close all, clc

% find input and output dirs
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir)
cd ..
cd outputs/
out_dir      = [pwd, '/entireModel/'];
in_dir_slab  = [pwd, '/slab/structures/'];
in_dir_UC    = [pwd, '/upCrust/structures/'];
out_dir_tl   = [pwd, '/tomolab/'];
cd ..
cd inputs/
in_dir_SR    = [pwd, '/stingray_structures/'];
cd(scriptDir)

addpath(out_dir); addpath(in_dir_slab); addpath(in_dir_UC); addpath(in_dir_SR);

%find the path
theUpCr      = which('upCrust_31-Mar-2025_UpCrust_Vp_12points_highSltzTerrane_AUG13_2024.xlsx.mat');   % 2D vp matrix for upper crust
theInt       = which('int_3D_31-Mar-2025_casieGP_bloch.mat'); % Interfaces
theStation   = which('srStation_29-Oct-2024_CG&SG_.mat');
theEvent     = which('srEvent_28-Oct-2024_CG&SG_.mat');
theElevation = which('srElevation_feb2_CGSG.mat');
theGeometry  = which('OR2012_srGeometry.mat');
theControl   = which('srControl_casExp.mat');

tlPickDir='/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/tlPicker/PickFiles_CG_SG';

phaseIn      = [{'P'} {'Pg'} {'Pn'} {'PmP'}];


%name the output model
unq_model_name  = 'test'; % generally put the region the experiment

disp('Making srModel ...')
run("make_srModel.m")
clc
disp('srModel saved')
disp('Making tlPert ...')
run('make_tlPert.m')
clc
disp('srModel saved')
disp('tlPert saved')

disp('Making tlArrival ...')
run('make_tlArrival.m')
clc
disp('srModel saved')
disp('tlPert saved')
disp('tlArrival saved')