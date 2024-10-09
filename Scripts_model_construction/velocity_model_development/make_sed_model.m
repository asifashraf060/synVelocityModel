 %% Script for developing a simplified model from a dense/detailed srModel
% make a 3D model with no variability along y axis/longitudal axis
% asif, March 28, 2022


% load the entire model
clear, close all
load('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/New_models_output/srModel_Cas2021_vF.mat');
%load('/Users/asifashraf/Talapas/tlOutput/220726_235058/srModel_it7.mat')
%make the model 3D to 2D
new_model = squeeze(srModel.P.u(:,1,:));
new_xPos  = srModel.xg; %xPos = position in horizontal scale
                        %new_model = squeeze(srModel.P.u(:,1,:));
new_xPos  = srModel.xg; %xPos = position in horizontal scale
                        %convert it back lat lon later
new_yPos  = srModel.yg; %yPos = changes in latitudinal directionconvert it back lat lon later
new_yPos  = srModel.yg; %yPos = changes in latitudinal direction

model_name = ['sed_2Dmat_', date, '_upCr_v0.mat'];

%manually manipulate the Vp
vp_1D_tb_11points = table2array(readtable('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/UpCrust_Vp_11points.xlsx'));
logical = 1;

%%%%%%%%%%%%%%%%%%% INDEX %%%%%%%%%%%%%%%%%%%%%%%
%% 1. Make a coarse model with 22 points throughout the extent
%% 2. Make another model with even coarser intervals (11 points throughout the model)
%% 3. Make a denser model with coarsely interpolated points
%% 4. Make a new srModel.interface.Z from smooth srModel.interface.elevation
%% 5. Now, cut the Vp, up to basement Z
%% 6. Filter the sedimentary velocity to make it smoother
%% 7. Need to fill the interfaces for rest of the x positions in west
%% 8. assign crustal velocities
%% 9. add the mantle velocities
%% 10. fill Vp a little bit more in the east direction
%% 11. fill Vp in a little bit more in the west direction too
%% 12. Create srModel structure
%%%%%%%%%%%%%%%%%%% INDEX %%%%%%%%%%%%%%%%%%%%%%%


%% 1. Make a coarse model with 22 points throughout the extent
%first decide how many points/1-D profile you want throughout the model
%then look at the number of x stations you have in the model
%then choose sampling interval based on the # of points you want

%extracting some 1-D profiles throughout the entire model
for i = 1:22
    a = new_model((50*i),:);
    vp(:,i) = a; %what I am extracting is vp, not u
    b = new_xPos(50*i); %extracting only x positions
    x(i,:)  = b;
    %extract the basement interface
    c = srModel.interface(1).elevation((50*i),:);
    base(:,i) = c;
    %extract the moho interface
    d = srModel.interface(2). elevation((50*i),:);
    mh(:,i) = d;
end

%manually manipulate the Vp
%vp_1D_tb = table2array(readtable('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/22_vp_lines.xlsx'));
%vp = vp_1D_tb;


xPos     = (x(:,1))';
basement = base(1,:);
moho     = mh(1,:);
zPos     = (linspace(-1, -50, 251))';
[X, Z]        = meshgrid(xPos, zPos);
[x, value, z] = griddata(xPos, zPos, vp, X, Z);
[depth, velocity]   = contourf(x, value, z, [1:.2:8]);
colormap(jet)
cc = colorbar;
ylabel(cc, 'Vp');
title('model with 22 points')
hold on
plot(xPos, basement, '.b')
plot(xPos, moho, '.r')
hold off

finer_xPos = linspace(xPos(1), xPos(end), 1000);
finer_zPos = (linspace(zPos(1), zPos(end), 500))';
%make interpolation between points
finer_vp   = interp2(xPos, zPos, vp, finer_xPos, finer_zPos, 'spline');

%% 2. Make another model with even coarser intervals (11 points throughout the model)
%first decide how many points/1-D profile you want throughout the model
%then look at the number of x stations you have in the model
%then choose sampling interval based on the # of points you want

%extracting some 1-D profiles throughout the entire model
for k = 1:11
    a3 = new_model((102*k),:);%102 is the sampling interval
    vp3(:,k) = a3; %what I am extracting is vp, not u
    
    %extracting only x positions in coarser interval (i.e., 11 points throughout the station)
    b3 = new_xPos(102*k);
    x3(k,:)  = b3;
    
    %extract the basement interface
    c3 = srModel.interface(1).elevation((102*k),:);
    base3(:,k) = c3;
    
    %extract the moho interface
    d3 = srModel.interface(2).elevation((102*k),:);
    mh3(:,k) = d3;
    
    %extract the basement interface in Z, not elevation
    e3  = srModel.interface(1).Z((102*k),:);
    base3_Z(:,k) = e3;
    
    %extract the moho interface in Z, not elevation
    f3 = srModel.interface(1).Z((102*k),:);
    moho3_Z(:,k) = f3;   
end


if logical == 1
    vp3  = vp_1D_tb_11points;
end


%make the nodes in x and z and interfaces from whatever came out from the loop
xPos3     = (x3(:,1))';
zPos3     = (linspace(-1, -50, 251))';
basement3 = base3(1,:);
moho3     = mh3(1,:);
basement3_z = base3_Z(1,:);
moho3_z     = moho3_Z(1,:);
            % in Z all I am getting are NaNs
            % need to develop my own Z later

[x11_grd, z11_grd] = meshgrid(xPos3, zPos3);

figure(99118), clf
[C, h] = contourf(x11_grd, z11_grd, (vp3), [0:.5:8.1]);
clabel(C, h)
colormap(jet)
colorbar
hold on
plot(-3.55, 0, 'vk')

%make nodes in x and z in finer intervals
%to get an interpolated Vp in finer scale, we will do scatteredInterpolant

%% 3. Make a denser model with coarsely interpolated points (11 points)
%now, make a function of scatteredInterpolant for a better interpolation of 11 points
xPos3_clV_sm = xPos3'; %clV means column vector and 'sm' means small
%create a column vector A for inputing all the longitudes
A = NaN(length(zPos3)*length(xPos3_clV_sm));
A = A(:,1);
  %making column vectors for the scatteredInterpolant function          
    A(1:length(zPos3)) = xPos3_clV_sm(1);
    A(1*length(zPos3)+1:2*length(zPos3)) = xPos3_clV_sm(2);
    A(2*length(zPos3)+1:3*length(zPos3)) = xPos3_clV_sm(3);
    A(3*length(zPos3)+1:4*length(zPos3)) = xPos3_clV_sm(4);
    A(4*length(zPos3)+1:5*length(zPos3)) = xPos3_clV_sm(5);
    A(5*length(zPos3)+1:6*length(zPos3)) = xPos3_clV_sm(6);
    A(6*length(zPos3)+1:7*length(zPos3)) = xPos3_clV_sm(7);
    A(7*length(zPos3)+1:8*length(zPos3)) = xPos3_clV_sm(8);
    A(8*length(zPos3)+1:9*length(zPos3)) = xPos3_clV_sm(9);
    A(9*length(zPos3)+1:10*length(zPos3)) = xPos3_clV_sm(10);
    A(10*length(zPos3)+1:11*length(zPos3)) = xPos3_clV_sm(11);
%assign a new variable to the column vector
xPos3_clV = A;
%now, do this thing for dpeth column
B = repmat(zPos3, length(xPos3), 1);
zPos3_clV = B;
%now, do this thing for vp3 matrix
C = NaN(length(zPos3)*length(xPos3_clV_sm));
C = C(:,1);
  %making column vectors for the scatteredInterpolant function
  %loop didn't work
    C(1:length(zPos3)) = vp3(:,1);
    C(1*length(zPos3)+1:2*length(zPos3)) = vp3(:,2);
    C(2*length(zPos3)+1:3*length(zPos3)) = vp3(:,3);
    C(3*length(zPos3)+1:4*length(zPos3)) = vp3(:,4);
    C(4*length(zPos3)+1:5*length(zPos3)) = vp3(:,5);
    C(5*length(zPos3)+1:6*length(zPos3)) = vp3(:,6);
    C(6*length(zPos3)+1:7*length(zPos3)) = vp3(:,7);
    C(7*length(zPos3)+1:8*length(zPos3)) = vp3(:,8);
    C(8*length(zPos3)+1:9*length(zPos3)) = vp3(:,9);
    C(9*length(zPos3)+1:10*length(zPos3)) = vp3(:,10);
    C(10*length(zPos3)+1:11*length(zPos3)) = vp3(:,11);    
%assign a new variable to the column vector
vp3_clV = C;
%make a function of scatteredInterpolant
F = scatteredInterpolant(xPos3_clV, zPos3_clV, vp3_clV, 'linear');
%now, make a finer interpolated vp3
%to do that, we need to make finer x and z nodes
finer_xPos3  = (linspace(xPos3(1), xPos3(end), 1000))';
finer_zPos3  = flip((linspace(-50, 0, length(finer_xPos3)))');
[x_fn, y_fn] = meshgrid(finer_xPos3, finer_zPos3);
%apply the scattered interpolant function to make a finer_vp3
finer_vp3 = F(x_fn, y_fn);

%% 4. Make a new srModel.interface.Z from srModel.interface.elevation
%first, make new interfaces for finer_xPos (i.e., finely grided x stations from coarsely interpolated model)
%basement
basement_interp   = interp1(xPos3, basement3, finer_xPos, 'cubic');
%do the same thing for moho
moho_interp       = interp1(xPos3, moho3, finer_xPos, 'cubic');
%interpolated data has NaN, because the loop in the previous section is not starting from first point
%so what we are doing here is replacing the NaN with the first value of the
%basement/moho
basement_interp(isnan(basement_interp)) = basement_interp(51); %51 is the first interpolated point from the loop
moho_interp(isnan(moho_interp))         = moho_interp(51);

%so, you have the basement and moho at the correct interval
%now, its time interpolate/extract elevation
%to bring the interface from elevation field to Z field

            %but, you are lacking the array of y positions
            %to apply griddata function on elevation, you need it
            %so, building it ---
            yPos_4oneLine = [zeros(length(finer_xPos3))+mean(new_yPos)];
            yPos_4oneLine = (yPos_4oneLine(:,1));
            %make a meshgrid to apply griddata
            [grd_x, grd_y] = meshgrid(finer_xPos3, yPos_4oneLine);
%now,
elevation_griddata = griddata(srModel.interface(1).X, srModel.interface(1).Y, ...
                                     srModel.elevation, grd_x, grd_y);
% make a nice array of elevation data
elevation_data     = elevation_griddata(1,:);
%FINALLY, its time to create your own Z field for elevation
basement_Z = basement_interp - elevation_data;
moho_Z     = moho_interp - elevation_data;
%check if it's correct
% figure(11), clf
% plot(finer_xPos, basement_interp, '.-b')
% hold on
% plot(finer_xPos, moho_interp, '.-r')
% plot(finer_xPos, elevation_data, '.-k')
% plot(finer_xPos, basement_Z, '-c')
% plot(finer_xPos, moho_Z, '-m')
% legend('basement_E', 'moho_E', 'elevation', 'basement_Z', 'moho_Z', 'Location', 'SouthWest', 'FontSize', 12)
% title('basement and moho made from only 11 points of finer model', 'FontSize', 16)
% xlabel('X');
% ylabel('depth (m)');


% Make a 3-D srModel.interface.Z which will account for elevation at every line along y axis
% for i_3D_y = 1:length(new_yPos)
%     
%     yPos_loop = [zeros(length(finer_xPos3))+(new_yPos(i_3D_y))];
%     yPos_loop = (yPos_loop(:,1));
%     [grd_x_loop, grd_y_loop] = meshgrid(finer_xPos3, yPos_loop);
%     elevation_griddata_loop = griddata(srModel.interface(1).X, srModel.interface(1).Y, ...
%                                      srModel.elevation, grd_x_loop, grd_y_loop);
%     elevation_data_loop     = elevation_griddata_loop(1,:);
%     
%     basement_Z_loop         = basement_interp - elevation_data_loop;
%     moho_Z_loop             = moho_interp - elevation_data_loop;
% 
%     basementZ_3D(:,i_3D_y)  = basement_Z_loop;
%     mohoZ_3D(:,i_3D_y)      = moho_Z_loop;
% end
% 
% figure(99123), clf
% [x_meshgrid, y_meshgrid] = meshgrid(new_yPos, finer_xPos3);
% plot3(x_meshgrid, y_meshgrid, basementZ_3D)
% hold on
% plot3(x_meshgrid, y_meshgrid, mohoZ_3D, 'y')
% grid on
% 
% figure(99124), clf
% [X, Y] = contourf(x_meshgrid, y_meshgrid, basementZ_3D, [-35:1:-5]);
% colorbar


%% 5. Now, cut the Vp, up to basement Z
%now make a matrix filling up rest of the position of xPos
fill_xPos3 = -144.2:.2:-133.8;
fill_vp3   = repmat(finer_vp3(:,1), [1, length(fill_xPos3)]);
%now concat
final_vp3  =horzcat(fill_vp3, finer_vp3);
final_xPos3 = horzcat(fill_xPos3, finer_xPos);
%now, its time to extract another form of basement and moho
%to do that, create meshgrid
[finer_x, finer_z] = meshgrid(finer_xPos, finer_zPos3);
[final_x, final_z] = meshgrid(final_xPos3, finer_zPos3);
basement_Z_repMat  = repmat(basement_Z, [length(basement_Z), 1]);
%interpolate basement in same interval
final_basement_Z  = griddata(finer_x, finer_z, basement_Z_repMat, final_x, final_z);
final_basement_Z  = final_basement_Z(1,:);
%do the same thing for moho
moho_Z_repMat     = repmat(moho_Z, [length(moho_Z), 1]);
final_moho_Z      = griddata(finer_x, finer_z, moho_Z_repMat, final_x, final_z);
final_moho_Z      = final_moho_Z(1,:);

%finally, make a loop to cut it
for ijkl = 1:length(final_z)
    df = find(final_z(:,ijkl)>final_basement_Z(ijkl) == 1);
    fg = final_vp3(:,ijkl);
    gh = fg(df); % cut vp upto the basement
    hj = NaN((1000-length(gh)),1); %fill the rest of the column with NaN
    jk = vertcat(gh, hj);
    sed_vp(:,ijkl) = jk; %sed_vp = Vp of layers on top of slab
end

%% 3D fill the basement and Moho for the rest of the x positions in west

% for i_3D_fill = 1:length(new_yPos)
% yPos_loop2 = [zeros(length(fill_xPos3))+new_yPos(i_3D_fill)];
% yPos_loop2 = yPos_loop2(:,1);
% [grd_x_loop2, grd_y_loop2] = meshgrid(fill_xPos3', yPos_loop2);
% 
% qaz = basementZ_3D(:,i_3D_fill);
% basementZ_3D_fill = repmat(qaz(1), [1, length(fill_xPos3)]);
% basementZ_3D_fill_fill(:,i_3D_fill) = basementZ_3D_fill;
% 
% waz = mohoZ_3D(:,i_3D_fill);
% mohoZ_3D_fill = repmat(waz(1), [1, length(fill_xPos3)]);
% mohoZ_3D_fill_fill(:, i_3D_fill) = mohoZ_3D_fill;
% end
% 
% %add the 3D interfaces together
% basementZ_3D_final = vertcat(basementZ_3D_fill_fill, basementZ_3D);
% mohoZ_3D_final = vertcat(mohoZ_3D_fill_fill, mohoZ_3D);

%% 6. Filter the sedimentary velocity to make it smoother

%now, fill all the NaNs of sed_vp with the deepest value
for ijklm = 1:length(sed_vp)
    mko   = sed_vp(:,ijklm);
    mko(find(isnan(mko))) = [];
    bhu   = sed_vp(:,ijklm);
    bhu(find(isnan(bhu))) = mko(end);
    sed_vp_woNaNs(:,ijklm) = bhu;
end


%do the gauss filtering for the sedimentary velocity
sed_vp_gfilt = imgaussfilt(sed_vp_woNaNs);
sed_vp_gfilt50= imgaussfilt(sed_vp_woNaNs, 50);
sed_vp_gfilt25= imgaussfilt(sed_vp_woNaNs, 25);
sed_vp_gfilt10= imgaussfilt(sed_vp_woNaNs, 10);
sed_vp_gfilt5= imgaussfilt(sed_vp_woNaNs, 5);


figure(2), clf
[C,h] = contourf(flipud(sed_vp_gfilt5));
colormap(jet)
clabel(C, h)

figure(3), clf
[C,h] = contourf(final_x, final_z, (sed_vp_gfilt25),[0:.2:8.1]);
colormap(jet)
clabel(C, h)
hold on
plot(-3.55, 1.5, 'vk')

sed_2Dmat.vp = sed_vp_gfilt25;
sed_2Dmat.xPos = final_xPos3;
sed_2Dmat.zPos = finer_zPos3;

save(['/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/', model_name],'sed_2Dmat')