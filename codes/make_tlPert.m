% Script to make tlPert
% Asif, Apr 21

%% Load srStation and srEvent

load(theModel)
load(theGeometry)
srStation=load_srStation(theStation, srGeometry);
srEvent=load_srEvent(theEvent, srGeometry);

%% Make tlPert

% MODEL parameteres
Model_minz   = min(srModel.zg);  %km
% PERT parameters
Pert_dx = 2;       %km
Pert_dy = 2;       %km
Pert_dz = 2;       %km
Pert_unc = 0.05;   %km
Pert_moho_dx = 2;       %km
Pert_moho_dy = 2;       %km
Pert_moho_unc = 0.05;   %km
tf_Pert_anis = 0;       %km

% make tlPert structure
tlPert.U.P.x = srModel.xg(1):Pert_dx:srModel.xg(end);
tlPert.U.P.y = srModel.yg(1):Pert_dy:srModel.yg(end);
tlPert.U.P.z = [srModel.zg(1):-Pert_dz:srModel.zg(end)]';  % Careful z is negative

% tlPert must have a node at the ends of the models.
tlPert.U.P.x = [tlPert.U.P.x(1:end-1) srModel.xg(end)];
tlPert.U.P.y = [tlPert.U.P.y(1:end-1) srModel.yg(end)];
tlPert.U.P.z = [ tlPert.U.P.z(1:end-1); Model_minz ];
tlPert.U.P.nx = length(tlPert.U.P.x);
tlPert.U.P.ny = length(tlPert.U.P.y);
tlPert.U.P.nz = length(tlPert.U.P.z);

tlPert.U.P.unc = Pert_unc*ones(tlPert.U.P.nx,tlPert.U.P.ny,tlPert.U.P.nz);

tlPert.I.id = srModel.interface(1).id;
tlPert.I.name = srModel.interface(1).name;
tlPert.I.x = tlPert.U.P.x(1):Pert_moho_dx:tlPert.U.P.x(end);
tlPert.I.y = tlPert.U.P.y(1):Pert_moho_dy:tlPert.U.P.y(end);

tlPert.I.x = [tlPert.I.x(1:end-1) srModel.xg(end)];
tlPert.I.y = [tlPert.I.y(1:end-1) srModel.yg(end)];
nx = length(tlPert.I.x);
ny = length(tlPert.I.y);

tlPert.I.unc = Pert_moho_unc*ones(nx,ny);
tlPert.I.unc = tlPert.I.unc;

thePert = [out_dir_tl, 'tlPert_', date, '_', unq_model_name, '.mat'];
save(thePert,'tlPert')
tlPert = load_tlPert(thePert, srModel);
save(thePert,'tlPert')