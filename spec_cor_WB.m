function [datpbin_a_corr,datpbin_c_corr]=spec_cor_WB(datpbin,band,cols,cbreak,abreak)   
    %%Spectral discontinuity correction
    % seperate to get beam.c & beam.a
    bandsep=(cols/2); %find the break between beam.c and beam.a
    datpbin_c=datpbin(:,2:bandsep);
    datpbin_a=datpbin(:,bandsep+1:end-1); %last column is all 0's
    mdate=datpbin(:,1);

    %cbreak & abreak are the two wavelenths where the spectra breaks for
    %beam_c and beam_a respectively 
    %


    %for beam.c
    %compute slope AFTER the break in spectra
    slope=(datpbin_c(:,cbreak(2)+1)-datpbin_c(:,cbreak(2)))/(band(cbreak(2)+1,1)-band(cbreak(2),1)); %band 1st column for beam.c no.
    %compute where the data should be
    intercept=datpbin_c(:,cbreak(2))-slope*band(cbreak(2),1);
    %compute the difference btwn where
    per_c=slope.*band(cbreak(1),1)+intercept;
    diff_c=per_c-datpbin_c(:,cbreak(1));
    datpbin_c_corr=datpbin_c;
    datpbin_c_corr(:,1:cbreak(1))=datpbin_c(:,1:cbreak(1))+diff_c*ones(1,cbreak(1));
    
    %for beam.a
    slope=(datpbin_a(:,abreak(2)+1)-datpbin_a(:,abreak(2)))/(band(abreak(2)+1,2)-band(abreak(2),2)); %band 2nd column for beam.a band no.
    intercept=datpbin_a(:,abreak(2))-slope*band(abreak(2),2);
    per_a=slope.*band(abreak(1),2)+intercept;
    diff_a=per_a-datpbin_a(:,abreak(1));
    datpbin_a_corr=datpbin_a;
    datpbin_a_corr(:,1:abreak(1))=datpbin_a(:,1:abreak(1))+diff_a*ones(1,abreak(1));


return