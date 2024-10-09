%to cut only the upper crustal portion of any inverted srModel
close all, clear all

%% INPUT
%load the inverted srModel
load('/Users/asifashraf/Talapas/tlOutput/230727_000634/srModel_it7.mat')
%showing the model
figure(1), clf
imagesc(squeeze((1./srModel.P.u (:,round(abs(length(srModel.yg)/2)) ,:)))');
title('loaded inverted srModel')
%input an interface structure
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/output_int_and_UpCr_struc/int_3D_28-Jul-2023int3D_26-Jul-2023int3D_pd16to18_ps01B_td1617_mcsDpCnv_thckTstd_DOWN_3to30.mat');
int = int_3Dmat;
%assign constant velocities
top_crustal_vp    = 6.1;
bottom_crustal_vp = 6.8;
mantle_vp         = 8.1;
max_depth         = 90;

%OUTPUT
out_dir = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/output_new_models';
model_name = append('srModel_', date, '_frmINV_ppgpnRun_down3');
disp('files are loaded')

%CALCULATION
%make the interface same as the model extent
[int_X, int_Y]     = meshgrid(int.yPos, int.xPos);
[model_X, model_Y] = meshgrid(srModel.yg, srModel.xg);
basement_Z         = interp2(int_X, int_Y, int.basement', model_X, model_Y);
moho_Z             = interp2(int_X, int_Y, int.moho', model_X, model_Y);

%Loop to assign crustal and moho velocity
%to use the inverted model's own interface
disp('Cutting inverted model velocities...')
moho     = moho_Z;
basement = moho+6;
if min(min(moho))<(max_depth*(-1))
    error('Maximum depth needs to be increased')
end
diff     = abs(srModel.zg(1) - srModel.zg(2));
zg_add   = flip(((max_depth*(-1)):diff:(min(srModel.zg) - diff))');
zg       = vertcat(srModel.zg, zg_add);
for i = 1:length(srModel.yg)
    up2D      = squeeze(srModel.P.u (:,i,:));
    vp2D      = 1./up2D;
    int_1_2D  = squeeze(basement (:,i,:));
    int_2_2D  = squeeze(moho (:,i,:));
    for K = 1:length(srModel.xg)
        vp1D           = vp2D(K,:);
        base_start     = find(zg<int_1_2D(K));
        moho_id        = find(zg<int_2_2D(K));
        vp1D(moho_id)  = mantle_vp;
        base_id        = [(base_start(1)):1:(moho_id(1))];
        basement_vp    = (linspace(top_crustal_vp, bottom_crustal_vp, length(base_id)));
        vp1D(base_id)  = basement_vp;
        vp_2D(K,:)     = vp1D;
    end
    vp_3D_model(:,i,:) = vp_2D;
end

figure(2), clf
[model_xz, model_zx]   = meshgrid(zg, srModel.xg);
meanY = round(abs(length(srModel.yg)/2));
abc = (squeeze(vp_3D_model(:,meanY,:)))';
contourf(model_xz, model_zx, abc' ,[0:.2:8.1])
colormap(jet)
title('model made from imported model interfaces', 'FontSize', 16)

%% Change existing srStructure

%first, remove previous interface structures
srModel = rmfield(srModel, 'interface');
%second, assign new values to those structures
srModel.P.u      = 1./vp_3D_model;
srModel.zg       = zg;
srModel.nz       = length(srModel.zg);
minz             = min(srModel.zg);
srModel.ghead(5) = length(srModel.zg);

moho_interp      = (interp2(int_3Dmat.xPos, int_3Dmat.yPos', int_3Dmat.moho,...
                            srModel.xg, srModel.yg'))'; 
moho_interp_elv  = (interp2(int_3Dmat.xPos, int_3Dmat.yPos', int_3Dmat.moho_elv,...
                            srModel.xg, srModel.yg'))';

srModel.interface(1).id                         = 1;
srModel.interface(1).name                       = {'moho'};
[srModel.interface(1).Y,srModel.interface(1).X] = meshgrid(srModel.yg, srModel.xg);
srModel.interface(1).Z                          = moho_interp;
srModel.interface(1).elevation                  = moho_interp_elv;
%fixing the elevation
elevation = srModel.elevation;
save((append(out_dir, '/', model_name, '.mat')), 'srModel')

theElevation = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srElevation_s3k2k_grid123.mat';
theControl   = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/srControl_casExp.mat';
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Stingray_structures/OR2012_srGeometry.mat')
%load the structures through stingray
disp('Saving srModel...')
theModel     = (append(out_dir, '/', model_name, '.mat'));
srElevation  = load_srElevation(theElevation, srGeometry);
srControl    = load_srControl(theControl);
srModel      = load_srModel(theModel, srControl, srGeometry, srElevation);

figure(33), clf
plot3(srModel.interface(1).X, srModel.interface(1).Y, srModel.interface(1).Z, '-y')
hold on
plot3(srModel.interface(1).X, srModel.interface(1).Y, srModel.interface(1).elevation, '-r')
grid on

figure(34), clf
[ZZ, XX] = meshgrid(srModel.xg, srModel.zg);
vp = 1./srModel.P.u;
meanY = round(abs(length(srModel.yg)/2));
xSection = (squeeze(vp(:,meanY,:)))';
xInt_Z   = squeeze(srModel.interface(1).Z(:,meanY));
xInt_elv = squeeze(srModel.interface(1).elevation(:,meanY));
contourf(ZZ, XX, xSection)
hold on
plot(srModel.xg, xInt_Z, '.r')
plot(srModel.xg, xInt_elv, '.b')

srModel.elevation = elevation;
disp(srModel)
save(theModel, 'srModel')
disp('New srModel saved')

