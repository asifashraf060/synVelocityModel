clear all, close all
tb_TOC = table2array(readtable('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/interface_text_files/Carbotte-TOC-twtt-PD16_17_18-PS01B_TD1617_LL.txt'));
tb_SF  = table2array(readtable('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/interface_text_files/Carbotte-SF-twtt-final_LL.txt'));


%% TOC & SF time

ln_toc = tb_TOC(:,1); lt_toc = tb_TOC(:,2); tm_toc = tb_TOC(:,3);
ln_sf  = tb_SF(:,1);  lt_sf  = tb_SF(:,2);  tm_sf  = tb_SF(:,3);

figure(1),clf
plot3(ln_toc, lt_toc, tm_toc, '.r')
grid on
hold on
plot3(ln_sf, lt_sf, tm_sf, '.b')
set(gca, 'Zdir', 'reverse')


% interpolate


tm_sf_int = griddata(ln_sf, lt_sf, tm_sf, ln_toc, lt_toc, 'nearest');

if ~isempty(find(isnan(tm_sf_int)==1))
   warning(append('Interpolating from seafloor interface contains ', string(length(find(isnan(tm_sf_int)==1))) ,' nans')) 
else
   disp('There is no nan in the interpolation from seafloor interface')
end


%element 

only_toc_tm = (tm_toc./2) - (tm_sf_int./2);

only_toc_tm(find(only_toc_tm<0)) = 0;


figure(2),clf
plot3(ln_toc, lt_toc, tm_toc./2, '.r')
grid on
hold on
plot3(ln_toc, lt_toc, tm_sf_int./2, '.b')
plot3(ln_toc, lt_toc, only_toc_tm, '.g')
set(gca, 'Zdir', 'reverse')

fid = fopen(append(append('/Users/asifashraf/Documents/MATLAB/Cascadia Analysis/Initial_Model/InitialModel_development/input_files/', 'time_TOC_SFsubtracted.txt')), 'w');

for i = 1:length(only_toc_tm)
    fprintf(fid, '%f %f %f\n', ln_toc(i), lt_toc(i), only_toc_tm(i));
end

fclose(fid);