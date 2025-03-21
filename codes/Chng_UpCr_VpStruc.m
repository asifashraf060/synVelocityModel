%original excel with 11 1-D Vp points
vp_ar = table2array(readtable('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/UpCrust_Vp_11points.xlsx'));
%Range of column for Sltz Vp
R = (8:11);
%Path to newly developed Vp structure
vp_v11_path = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/UpCrust_vp_V11.xlsx';
vp_v22_path = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/UpCrust_vp_V22.xlsx';
vp_v33_path = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/UpCrust_vp_V33.xlsx';
vp_v44_path = '/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/BIG_matrix/UpCrust_vp_V44.xlsx';

%choose the column for only Siletz Vp model
sltz_vp = vp_ar(:,R);

%V11 -- +1 than the Vp
v11_add  = sltz_vp + 1;

%V22 -- -1 than the Vp
v22_add  = sltz_vp - 1;

%V33 -- +0.5 in upper 20%
A                        = size(sltz_vp);
up_20p                   = A(1) * 0.2;
up_20p                   = round(up_20p);
sltz_vp_up20p            = sltz_vp((1:up_20p),:);
sltz_vp_up20p_v33        = sltz_vp_up20p + 0.5;
sltz_vp_md               = sltz_vp;
sltz_vp_md((1:up_20p),:) = sltz_vp_up20p_v33;
v33_add                  = sltz_vp_md;

%v44 -- +0.5 in lower 80%
sltz_vp_bt80p                  = sltz_vp(((up_20p+1):end),:);
sltz_vp_bt80p_v33              = sltz_vp_bt80p + .5;
sltz_vp_md                     = sltz_vp;
sltz_vp_md(((up_20p+1):end),:) = sltz_vp_bt80p_v33;
v44_add                        = sltz_vp_md;

%make full array with modified vp columns
vp_md = vp_ar;
vp_md(:,R) = v11_add;
v11 = vp_md;
vp_md = vp_ar;
vp_md(:,R) = v22_add;
v22 = vp_md;
vp_md = vp_ar;
vp_md(:,R) = v33_add;
v33 = vp_md;
vp_md = vp_ar;
vp_md(:,R) = v44_add;
v44 = vp_md;

writematrix(v11, vp_v11_path, 'Sheet', 1)
writematrix(v22, vp_v22_path, 'Sheet', 1)
writematrix(v33, vp_v33_path, 'Sheet', 1)
writematrix(v44, vp_v44_path, 'Sheet', 1)