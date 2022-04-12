close all
clear all
clc

%% ACS processing -  Step 3
% 24/4/2021
% Ben Lowin

% This part of the script corrects for spectral discontinuity
% Corrects for the effect of tempurature
% Calculates Chlorophyll


%% User inputs
% You need to fill out this section with the relavant information. 
%

%make sure you are in the right directory

%load in the acs_cruiseID_binned.mat
load acs_SKQ_2021_jun_binned.mat


%% Spectral discontinuity correction
% First column in acs data (datpbin) is time, columns 2-83 are beam
% attenuation for the 82 unique wavelengths, % columns 84-165 are
% absorption for the 82 unique wavelengths and the last column is blank,
% zeros only. 
%
% The program will randomly choose a wavelenght to display, this can be
% re-run to find other wavelenghts.
% there are three section that should be run when doing final processing. 
% they care commented out currently. They are to look an change the break
% points.



find(~isnan(datpbin(:,159))); %to ensure you don't choose a NAN in the plots below
p=randsample(ans,1); %name the randomly chosen timepoint p

%Find where the spectral discontinuity is for the c-beam
%{
***COMMENTED OUT FOR AT_SEA PROCESSING, NOT NECESSARY
 figure
 for i=1:10
 plot(datpbin(p+i,2:cols/2),'o-')
 hold on
 end
%Zoom in on middle of spectra. Based on this plot, write down the break in the attenuation (c) spectra
%}

cbreak=[44 45]; %c-break point, check for final processing

%Find where the spectral discontinuity isfor the a-beam
%{
***COMMENTED OUT FOR AT_SEA PROCESSING, NOT NECESSARY
 figure
 for i=1:10
 plot(datpbin(p+i,(cols/2+1):end-1),'o-')
 hold on
 end
Zoom in on middle of spectra. Based on this plot, write down the break in the attenuation (c) spectra
%}
abreak=[45 46];  %a-break point, check for final processing


%This shifts the spectra vertically before and after a gap in wavelength space
[datpbin_a_corr,datpbin_c_corr]=spec_cor_WB(datpbin,band,cols,cbreak,abreak);
%datpbin_a_corr is now the absorption data with spectral discontinuity corrected
%datpbin_c_corr is now the absorption data with spectral discontinuity corrected

%check to see if it worked properly. We now plot spectra at random times.
figure
subplot(2,1,1)
plot(datpbin(p,2:cols/2),'o-') 
hold on; plot(datpbin_c_corr(p,:),'o-')
legend('uncorrected','corrected')
title('Spectral discontinuity correction for the c beam')

subplot(2,1,2)
plot(datpbin(p,(cols/2+1):end-1),'o-')
hold on; plot(datpbin_a_corr(p,:),'o-')
legend('uncorrected','corrected')
title('Spectral discontinuity correction for the a beam')

%% RE-SAMPLE A & C Bean
% changes from evey wavelent to evey other wavelengh
% check to make sure you have the correct number of wavebands

if length(datpbin_a_corr(1,:)) ~= length(band(:,1))
    error('wavelength mismatch. check bandname file')
end
band_new=[400:2:750]; %these are the new 2nm wavebands we are interpolating to
[dat_a,dat_c]=center_band_WB(datpbin_a_corr,datpbin_c_corr,band,band_new);


%% Residual Temperature correction
% Subtract a wavelength-dependent temp-salinity signal based on effects of pure water (Sullivan table)
% Read in the Sullivan Table to get Ds coefficient
Fai_S=textread('Sullivan_table.txt','',-1,'headerlines',1); %6th is what we want
[ap_cor,cp_cor]=resT_corr_WB(dat_a,dat_c,Fai_S); %

max_band=round(min(max(band))); %find the max wavebands for a & c, take the lower one, and round it
max_interp=max(find(band_new<=max_band)); %find bands that we actually have data for
ap_cor=ap_cor(:,1:max_interp); %remove columns where we have no data
cp_cor=cp_cor(:,1:max_interp); %remove columns where we have no data

%% Set up the save out structure
%
%

acs.ap=ap_cor; %particulate absorption data for all wavelenghts
acs.cp=cp_cor; %particulate attenuation data for all wavelenghts
acs.mdate=datpbin(:,1); %time (mdate)
acs.jd=acs.mdate-ref_date; %time (jd)
acs.wavelength=band_new(1:max_interp); %acs wavelenghts


%% Calculate Chl (using line height technique) 
% pulls ot the corosponding wavelenght data srounding the peak made by
% chlrophyll at 676nm. 
%

a650=acs.ap(:,find(acs.wavelength==650));
a676=acs.ap(:,find(acs.wavelength==676));
a714=acs.ap(:,find(acs.wavelength==714)); %this is 714nm!!!!
a716=acs.ap(:,find(acs.wavelength==716)); %this is 716nm!!!!

%avrages 714 and 716 to find 715 nm.
for i=1:length(acs.ap) 
    a715(i)=nanmean([a714(i),a716(i)]);
end

a715=a715'; % inverts to have in correct oreintation
abl=((a715-a650)/(715-650)).*(676-650)+a650;

acs.lha676=a676-abl; %adds to the acs structure - line height absorbstion
acs.chl=acs.lha676/0.0126; %coefficient from Chl paper ??? Which one? (Gaff et al. 2016)?? Guess?

%plot the final Chl dataset
figure
plot(acs.mdate-ref_date,acs.chl,'.')
title('Chl timeseries [mg m-3]')
xlabel('julian day')
ylabel('Chl [mg m-3]')

save(horzcat('acs_',cruiseID,'_processed','.mat'),'acs','ref_date');

clearvars a650 a676 a715 a714 a716
