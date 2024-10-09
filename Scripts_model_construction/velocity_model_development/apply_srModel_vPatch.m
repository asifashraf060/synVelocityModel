function [matrix] = apply_srModel_vPatch(srModel, vp_input)

    
    [xz, zx]    = meshgrid(srModel.xg, srModel.zg);
    vp2D        = squeeze(1./srModel.P.u (:,round(length(srModel.yg)/2),:));
    basement_2D = squeeze(srModel.interface(1).Z (:,round(length(srModel.yg/2)),:));
    
    %plot the middle yline
    figure(110), clf
    contourf(xz, zx, vp2D', [min(min(vp2D)):.2:max(max(vp2D))])
    colormap('jet')
    colorbar
    title('Vp Patch: select 2 points above the slab', 'FontSize', 16);
    [px, pz]= ginput(2);
    
    %indices for all 1-D profiles
    ind1              = find(srModel.xg>px(1) & srModel.xg<px(2));
    %get all the x points
    xP                = srModel.xg(ind1);
    basement_selected = basement_2D(ind1)+8;
    
    figure(220), clf
    plot(px, pz, '-or', 'DisplayName', 'Selected Points')
    hold on
    plot(xP, basement_selected, '-b', 'DisplayName', 'Basement')
    ylabel('Depth (km)')
    xlabel('xg')
    legend('Location', 'southwest')
    
    %constant thickness
    d = abs(basement_selected(1) - pz(1));
    
    %loop to change the vp in upward direction of the basement 
    for i = 1:length(ind1)
        id   = ind1(i);
        vp1D = vp2D(id,:);
        base = basement_selected(i);
        top  = base + d;
        mid  = find(srModel.zg>base & srModel.zg<top);
        vp1D(mid)  = vp_input;
        vp2D(id,:) = vp1D;
    end
    
    figure(330), clf
    contourf(xz, zx, vp2D', [min(min(vp2D)):.2:max(max(vp2D))])
    colormap('jet')
    colorbar
    title('Changed Vp model', 'FontSize', 16);
    
    vp3D = permute(repmat(vp2D, 1, 1, length(srModel.yg)), [1,3,2]);
    matrix = 1./vp3D;

end