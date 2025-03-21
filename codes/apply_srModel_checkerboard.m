function [matrix] = apply_srModel_checkerboard(srModel, prcntDv, waveDimension, condition, upCrustModel)

z1 = max(srModel.zg); z2 = min(srModel.zg);

w  = waveDimension;

v = 1./srModel.P.u;

basement_interp = srModel.interface(1).Z + 6;   moho_interp = srModel.interface(1).Z;

[~,~,k1] = xyz2ijk(0,0,z1,srModel.ghead); [~,~,k2] = xyz2ijk(0,0,z2,srModel.ghead);

xrange = length(v(:,1,1)); yrange = length(v(1,:,1)); zrange = length(v(1,1,:));


% sine wave with 1 to -1 value amp
sinx = sin(2*pi*(1:(xrange))/w(1)); % ( sin( 2*pi*(1:(xrange))/w(1) ) + 1 )/2; %( sin( 2*pi*(9:(xrange+8))/w(1) ) + 1 )/2;
siny = sin(2*pi*(6:(yrange+5))/w(2)); %( sin( 2*pi*(6:(yrange+5))/w(2) ) + 1 )/2; % ( sin( 2*pi*(18:(yrange+17))/w(2) ) + 1 )/2;
sinz = sin(2*pi*(1-z1+4:(zrange-z1+4))/w(3)); % ( sin( 2*pi*(1-z1+4:(zrange-z1+4))/w(3) ) + 1 )/2; % ( sin( 2*pi*(32:(zrange+31))/w(3) ) + 1 )/2;

for i = 1:xrange
    for j = 1:yrange
        for k = k1:k2
            v(i,j,k) = ( 1 - (prcntDv/3)*(sinx(i) + siny(j) + sinz(k))) * v(i,j,k);
        end
    end
end

if strcmp(condition, 'whole') == 1
    matrix = 1./v;
end

if strcmp(condition, 'onlyUpMantle') == 1
    for i = 1:length(srModel.yg)
        vp_2D_loop_uC = squeeze(upCrustModel(:, i ,:));
        vp_2D_loop_ck = squeeze(v(:,i,:));
        basement_loop = basement_interp(:,i);
        moho_loop     = moho_interp(:,i);
            for j = 1:length(srModel.xg)
                vp_1D_uC          = vp_2D_loop_uC(j,:);
                vp_1D_ck          = vp_2D_loop_ck(j,:);

                basement_loop2    = basement_loop(j);
                moho_loop2        = moho_loop(j);

                cut_index1        = find(srModel.zg<basement_loop2 & srModel.zg>moho_loop2);
                vp_crust          = linspace(6.1, 6.8, length(cut_index1));
                vp_1D_ck(cut_index1) = vp_crust;

                cut_index2        = 1:1:find(abs(srModel.zg - basement_loop2) == min(abs(srModel.zg - basement_loop2)));
                vp_1D_ck(cut_index2) = vp_1D_uC(cut_index2);

                vp_crust2D(j,:)   = vp_1D_ck;
            end
        vp_crust3D(:,i,:) = vp_crust2D;
    end
    matrix = 1./vp_crust3D;
end

if strcmp(condition, 'onlyUpCrust') == 1
    for m = 1:length(srModel.yg)
        vp_2D_loop    = squeeze(v(:,m,:));
        basement_loop = basement_interp(:,m);
        moho_loop     = moho_interp(:,m);
             for n = 1:length(srModel.xg)
                 vp_1D           = vp_2D_loop(n,:);
                 basement_loop2  = basement_loop(n);
                 moho_loop2      = moho_loop(n);
                 cut_index       = find(srModel.zg<moho_loop2);
                    if isempty(cut_index)
                        error('deepest moho is deeper than assigned highest moho depth')
                    end
                 vp_moho         = zeros(length(cut_index))+8.1;
                 vp_moho         = vp_moho(1,:);
                 vp_1D(cut_index)= vp_moho;
                 vp_moho2D(n,:)  = vp_1D;
             end
        vp_moho3D(:,m,:) = vp_moho2D;
    end
    matrix = 1./vp_moho3D;
end

end