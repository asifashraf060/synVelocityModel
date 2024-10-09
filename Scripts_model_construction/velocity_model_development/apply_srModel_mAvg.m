function srModel = srModel_mAvg(srModel, windowSize, direction)

    kernel = ones(1, windowSize) / windowSize;
    vp_3D  = 1./srModel.P.u;
    
    switch direction
        case 'UpDown'
    
            vp_3D_cnv = [];
             for i = 1:length(vp_3D)
    
                vp_2D = squeeze(vp_3D(i,:,:));
    
                vp_2D_cnv = [];
                for j = 1: length(vp_2D)
                    vp_1D = vp_2D(j,:);
                    vp_1D_cnv = conv(vp_1D, kernel, "same");
                    vp_2D_cnv(j,:) = vp_1D_cnv;
                end
              
                vp_3D_cnv(i,:,:) = vp_2D_cnv;
    
             end
    
        case 'LeftRight'
    
            vp_3D_cnv = [];
             for i = 1:width(vp_3D)
    
                vp_2D = squeeze(vp_3D(:,i,:));
    
                vp_2D_cnv = [];
                for j = 1: length(vp_2D)
                    vp_1D = vp_2D(j,:);
                    vp_1D_cnv = conv(vp_1D, kernel, "same");
                    vp_2D_cnv(j,:) = vp_1D_cnv;
                end
              
                vp_3D_cnv(:,i,:) = vp_2D_cnv;
    
             end
             
    end
    vp_3D_cnv_filt = imgaussfilt(vp_3D_cnv, 20);
    
    srModel.P.u = (1./vp_3D_cnv_filt);

end