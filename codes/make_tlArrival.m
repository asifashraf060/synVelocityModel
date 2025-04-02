% original--sarah
% modified--asif

%load srModel
load(theModel)
load(theGeometry)

%% 
tlPert=load_tlPert(thePert, srModel)
srStation=load_srStation(theStation, srGeometry)
srEvent=load_srEvent(theEvent, srGeometry)

%%
stationIn=[];
phaseIn_names = cell2mat(phaseIn);

%%
tlArrival = tlPick2tlArrival(srEvent, srStation, tlPickDir, stationIn, phaseIn); %stationIn, phaseOut, rlim, channelIn)

%% Correct error

tlArrival_name = ['tlArrival_', date,'_', phaseIn_names,'_' , unq_model_name, '.mat'];
save([out_dir_tl, tlArrival_name], 'tlArrival')
