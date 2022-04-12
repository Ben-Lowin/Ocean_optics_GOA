function [dat_a,dat_c]=center_band_WB(datpbin_a_corr,datpbin_c_corr,band,band_new)
    %%re-sample 2nm center band for beam.a&c
    dat_a=zeros(size(datpbin_a_corr,1),numel(band_new)); %make empty matrix
    
    dat_c=zeros(size(datpbin_a_corr,1),numel(band_new)); %make empty matrix
   
    for i=1:length(datpbin_a_corr)
        dat_a(i,:) = interp1(band(:,2),datpbin_a_corr(i,:)',band_new);
        dat_c(i,:) = interp1(band(:,1),datpbin_c_corr(i,:)',band_new);
    end
return