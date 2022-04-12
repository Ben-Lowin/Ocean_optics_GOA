close all
clear all
clc

%% ACS processing -  Step 2
% 14/4/2021
% Ben Lowin

% You will need to stop and fill out a spread sheet
    % to fill out the spread sheet use the 'ginput' comand
    %further instruction in Initial QC section - !read before you run!
%


%% User inputs
% You need to fill out this section with the relavant information.


%load in the acs_cruiseID_compile.mat
load acs_SKQ_2021_jun_compile.mat

% load the bad data spread-sheet
% should be names bad_data_CruiseID.xlsx
% cruiseID=ShipCode_year_MMM

bad_data_sheet=('bad_data_SKQ_2021_jun.xlsx');

% IF No bad data sheet yet - 
% change the 1 to a 0
Bad_data_made=1;

%% COMMENT OUT THE PREVIOUS ADD BLANKS
% GO do it
%

%% Match BB3 and ACS to identify blank periods
% first it rounds the time data, to the nearest second, so that they can be compared and removes
% inf. decimal. rounding currently goes towards the nearest intiger. 
% Identifies when blank and transtion were happening, from col 5 of the BB3
% data. 

%Round BB3 time into distinct seconds to match up with ACS time and
%removes the never ending decimal point. Rounds to the nearest second. 
bb3mdate=datbb3(:,1)+ref_date;
round(datevec(bb3mdate));
bb3mdate=datenum(ans);
%Round ACS time into distinct seconds to match up with BB3 time
round(datevec(mdate)); 
mdate=datenum(ans);


% finds the blank and transtion periods
blanks=find(datbb3(:,5)==1); %find when the 3-way valve was set to 'blank' mode
bb3blankmdate=bb3mdate(blanks); %dates when 3-way valve was set to 'blank' mode
trans=find(datbb3(:,5)==2); %find when the 3-way valve was transitioning between 'regular' and 'blank' mode
bb3transmdate=bb3mdate(trans);


% finds the common valuse between the time and time when the the blank was
% happening.
[blankmdate,blankind] = intersect(mdate,bb3blankmdate); %blankind are indices of ACS blankmdates
[transmdate,transind] = intersect(mdate,bb3transmdate); %transind are indices of ACS total-to-blank transition times

% plots a figure so that visual inspection of the blanks and transition
% periods can be done. 
f=figure;
f.Name='Visual inspection for blanks and transtions';
set(gcf,'OuterPosition',[100 400 1400 500]')
plot(mdate-ref_date,dat(:,159),'k.-'); %plotting a678
hold on
plot(mdate(blankind)-ref_date,dat(blankind,159),'r.') %plotting indices only should overlay blanks
plot(mdate(transind)-ref_date,dat(transind,159),'g.') %plotting indices only should overlay filter-to-blank transition
leg=legend('raw data','all blanks','all transitions');
title('Raw Absorption Data at 678nm')
ylabel('a678 [m-1]')
xlabel('julian day')
hold off

%differnt method, using consolidator, is not functional
%{
% new (potentially better) code??
% [x,y]=consolidator(mdate,dat(:,159));
% iblank=interp1(x,y,bbp.blankmdate);
% itrans=interp1(x,y,bbp.transmdate);
% hold on
% plot(bbp.blankmdate-ref_date,iblank,'r.')
% plot(bbp.transmdate-ref_date,itrans,'g.')
%[blankmdate,blankis,ib] = intersect(x,bbp.blankmdate); %ia1 is now the indices of ACS blanktimes
%[transmdate,transis,ib2] = intersect(x,bbp.transmdate); %ia2 is now the indices of ACS blanktimes
%}

%% Initial QC 
% Remove bad blanks, Remove obvious bad data, 
% This is done with the Figure made in the previous section.
% Add missing blanks - If blanks are missing, you need to find where they
% are and create new 'fake' blanks

%In accompanying excel spreadsheet "bad_data_cruiseID.xlsx", note down (in julien days) the start and end times for:
    %1)areas of obvious bad data. 
    %2)each bad blank. 
%use [x,y]=ginput for this. Makes it faster and more precise.

%SHIFT BLANK TIMES 
    %NOTE: If blanks need shifting slightly, use this code to do so
%{
% shift=datenum(0,0,0,0,0,60); %currently set to shift by 60 seconds
    %make sure to run both loops
    
% for i=1:length(bb3blankmdate)
%     bb3blankmdate(i)=bb3blankmdate(i)+shift;
% end

% for i=1:length(bb3transmdate)
%     bb3transmdate(i)=bb3transmdate(i)+shift;
% end

% bb3blankmdate=bb3mdate(blanks); %need this in case we run this code multiple times and overshift data
% bb3transmdate=bb3mdate(trans); %need this in case we run this code multiple times and overshift data
%}

    
%REMOVE OBVIOUS BAD DATA/BLANKS
%***NOTE: If a bad_data spreadsheet is not yet made, comment out the text
%below
if Bad_data_made==1
    [num,txt,~] = xlsread(bad_data_sheet);
    for i=1:length(num(:,1))
        baddata.start(i)=num(i,1);
        baddata.end(i)=num(i,2);
        baddata.reason(i,:)=txt(i+1,3)'; 
        inmdate=find(mdate-ref_date> baddata.start(i) & mdate-ref_date<baddata.end(i));
        mdate(inmdate)=nan; %remove from mdate
        dat(inmdate,:)=nan; %remove from entire acs dataset
    end
else
    error('Pause here and make bad_data.xlsx sheet')
end
    

%ADD MISSING BLANKS
%If Blanks are missing (e.g. if BB3 was turned off), use this code to add
%June 2020- Blanks
%{
addblank1=[datenum(2020,07,12,11,35,48):datenum(0,0,0,0,0,1):datenum(2020,07,12,11,38,15)]';    
addblank2=addblank1+datenum(0,0,0,0,17,20);
addblank3=addblank2+datenum(0,0,0,0,17,20);
addblank4=addblank3+datenum(0,0,0,0,17,20);
addblank5=addblank4+datenum(0,0,0,0,17,20);
addblank6=addblank5+datenum(0,0,0,0,17,20);
addblank7=addblank6+datenum(0,0,0,0,17,20);
addblank8=addblank7+datenum(0,0,0,0,17,20);
addblank9=addblank8+datenum(0,0,0,0,17,20);
addblank10=addblank9+datenum(0,0,0,0,17,20);
addblank11=addblank10+datenum(0,0,0,0,17,20);
addblank12=addblank11+datenum(0,0,0,0,17,20);
addblank13=addblank12+datenum(0,0,0,0,17,20);
%addblank14=addblank13+datenum(0,0,0,0,17,20);
addblank14=[datenum(2020,07,12,15,21,08):datenum(0,0,0,0,0,1):datenum(2020,07,12,15,28,04)]';
addblanksmdate=vertcat(addblank1,addblank2,addblank3,addblank4,addblank5,addblank6,addblank7,addblank8,addblank9,addblank10,addblank11,addblank12,addblank13,addblank14);
%add blanks are just dates
%}

%Blanks April 2021
%
%{
addblank1=[datenum(2021,04,24,20,28,19):datenum(0,0,0,0,0,1):datenum(2021,04,24,20,38,19)]';
addblank2=[datenum(2021,04,24,21,28,39):datenum(0,0,0,0,0,1):datenum(2021,04,24,21,38,39)]';
addblank3=[datenum(2021,04,26,17,57,50):datenum(0,0,0,0,0,1):datenum(2021,04,26,18,07,50)]';
addblank4=[datenum(2021,04,26,18,58,01):datenum(0,0,0,0,0,1):datenum(2021,04,26,19,08,01)]';
addblank5=[datenum(2021,04,26,19,58,13):datenum(0,0,0,0,0,1):datenum(2021,04,26,20,08,13)]';
addblank6=[datenum(2021,04,26,20,58,50):datenum(0,0,0,0,0,1):datenum(2021,04,26,21,08,50)]';
addblank7=[datenum(2021,04,29,01,40,39):datenum(0,0,0,0,0,1):datenum(2021,04,29,01,50,39)]';
addblank8=[datenum(2021,04,29,02,40,59):datenum(0,0,0,0,0,1):datenum(2021,04,29,02,50,59)]';
addblank9=[datenum(2021,04,29,03,41,36):datenum(0,0,0,0,0,1):datenum(2021,04,29,03,51,36)]';
addblank10=[datenum(2021,04,29,04,40,56):datenum(0,0,0,0,0,1):datenum(2021,04,29,04,50,56)]';
addblank11=[datenum(2021,04,29,05,41,25):datenum(0,0,0,0,0,1):datenum(2021,04,29,05,51,25)]';
addblank12=[datenum(2021,04,29,06,42,28):datenum(0,0,0,0,0,1):datenum(2021,04,29,06,52,28)]';
addblank13=[datenum(2021,04,26,21,58,53):datenum(0,0,0,0,0,1):datenum(2021,04,26,22,08,53)]';
addblank14=[datenum(2021,04,29,17,57,53):datenum(0,0,0,0,0,1):datenum(2021,04,29,18,07,53)]';
addblank15=[datenum(2021,04,29,18,58,59):datenum(0,0,0,0,0,1):datenum(2021,04,29,19,08,59)]';
addblank16=[datenum(2021,04,29,19,58,10):datenum(0,0,0,0,0,1):datenum(2021,04,29,20,08,10)]';
addblank17=[datenum(2021,04,29,20,59,08):datenum(0,0,0,0,0,1):datenum(2021,04,29,21,09,08)]';
addblank18=[datenum(2021,04,26,21,59,54):datenum(0,0,0,0,0,1):datenum(2021,04,26,22,09,54)]';


addblanksmdate=vertcat(addblank1,addblank2,addblank3,addblank4,addblank5,addblank6, addblank7,addblank8,addblank9 ...
    ,addblank10,addblank11,addblank12,addblank13,addblank14,addblank15,addblank16,addblank17);
%}

%un comment the following when doing blank adds
%{
round(datevec(addblanksmdate));
addblanksmdate=datenum(ans);

bb3blankmdate=vertcat(bb3blankmdate,addblanksmdate);
[t,i]=sort(bb3blankmdate);
bb3blankmdate=bb3blankmdate(i(:,1));
clearvars addblank1 addblank2 addblank3 addblank4 addblank5 addblank6 addblank7 addblank8 addblank9 addblank10 addblank11 addblank12 addblank13 addblank14
%}

%DO NOT COMMENT-OUT BELOW HERE
% Re-run the intersect in case bb3blankmdate & bb3transmdate are changed during QC steps above
[blankmdate,blankind] = intersect(mdate,bb3blankmdate); %blankind are indices of ACS blankmdates
[transmdate,transind] = intersect(mdate,bb3transmdate); %transind are indices of ACS total-to-blank transition times

%Plot the 'good' blanks. If any have been removed, they will still be red
ff=figure;
ff.Name='Visual inspection removing bad data';
set(gcf,'OuterPosition',[100 400 1400 500]')
plot(mdate-ref_date,dat(:,159),'k.-'); %plotting a678
hold on
plot(mdate(blankind)-ref_date,dat(blankind,159),'r.') %plotting indices only should overlay blanks
plot(mdate(transind)-ref_date,dat(transind,159),'g.') %plotting indices only should overlay filter-to-blank transition
plot(blankmdate-ref_date,dat(blankind,159),'b.') %plotting all of the 'good' data in blue
leg=legend('raw data','all blanks','all transitions','good blanks');
title('QC- Absorption Data at 678nm')
ylabel('a678 [m-1]')
xlabel('julian day')
hold off

%% Standardise data to blanks
% Find the blank values and interpolate between them. This interpolation is
% then subtracted from the observed.
% Total=particulate + water
% blank = water
% unfilltered = total
% rearangin this equation we can find just the particulate fraction for the bio-optics. 
% NOTE: acs_rem_blank_WB is not slecting a blank value for the last blank in
% the data set (July 2020). 
%
% This figure is very important as it checks that the blank interpolation is
% working correctly. An easy check for this is when going throught a front
% the blank value shoulc change. 
%
% This is then saved out so that we have the QC before bin



[datp,yi,x159,y159]=acs_rem_blank_WB(dat,mdate,blankmdate,ref_date,blankind,cols);
%datp is the matrix of particulate data (i.e. total - removed)
%yi is the interpolated blank values for all wavelengths (i.e. a_w and c_w)
%x159 and y159 are the blank values chosen at each blank period at ~a670nm 
%x159 and y159 are outputted so we can plot results and see how well the removal function worked

datblanks=dat(blankind,:);
yi=yi';


f3=figure;
f3.Name='Visual inspection intupilation of blanks';
set(gcf,'OuterPosition',[100 400 1400 500]')
plot(mdate-ref_date,dat(:,159),'k.-'); %plotting a678
hold on
plot(mdate(blankind)-ref_date,dat(blankind,159),'r.') %plotting indices only should overlay blanks
plot(mdate(transind)-ref_date,dat(transind,159),'g.') %plotting indices only should overlay filter-to-blank transition
plot(blankmdate-ref_date,dat(blankind,159),'b.') %plotting all of the 'good' data in blue
plot(x159,y159,'mo','markerfacecol','m'); % plotting the intepolation in magenta
plot(mdate-ref_date,yi(:,159),'m.'); 
leg=legend('raw data','all blanks','all transitions','good blanks','final blank values chosen','interpolation between chosen values');
title('Blank Interpolation - Absorption Data at 678nm')
ylabel('a678 [m-1]')
xlabel('julian day')
hold off
 
cd (resultsfolder)
save(horzcat('acs_',cruiseID,'_QC','.mat'), 'datp','ref_date','cols','cruiseID','acsfolder','band');
%% Remove blanks
% Using +/- 5 points we make sure to remove all of the blank related data. 
% currently this removes 10 seconds ahead and 20 seconds behind. 
% NOTE these do not have to be the same between the BB3 and the ACS
%           ACS will most likely be much shorter then the BB3 as it takes
%           much less time to flush.
% In this section you need to check that all of the blanks were identified and removed.  

for p=1:length(blankind)
    datp(blankind(p)-5:blankind(p)+5,:) = nan; %use 5 point threshold to make sure every filtered datapoint gets removed
end
for p=1:length(transind) %redo for all valve switchover times
    if transind(p)>length(datp(:,1))-180
       datp(transind(p)-30:transind(p)+1,:) = nan; %if its the last blank of the dataset, don't remove data after
    else    
    datp(transind(p)-10:transind(p)+20,:) = nan; %remove 10 seconds before every valve switchover, and 180 after
    end
end


f4=figure;
set(gcf,'OuterPosition',[100 350 1400 500]')
plot(datp(:,1),datp(:,159),'.','color',[0.7 0.7 0.7]);
hold on
line(dat(:,1),zeros(size(dat(:,1))))
title('Particulate Absorption at 678nm - Blanks removed')
xlabel('Julian day')
ylabel('a_{p}678 [m-1]')
hold on
plot(datp(:,1),datp(:,159),'k.');
legend('raw-filtered (ap676)','zero (blanks should run through here)','ap676 without blanks')
hold off

%% 1 MINUTE  BIN
% Bins the data a minute intervals using the consolidator function ( By
% John D'Errico)
% Look at this figure to see if the data has been correctly binned. Red
% points should plot on top of black. 

minute=1/24/60;
[x,y]=consolidator(mdate,datp(:,159),'nanmedian',minute);  %find out how many rows you'll need
datpbin=zeros(length(x),length(dat(1,:))); %initialize datpbin
for i=2:length(dat(1,:))-1
    [x,datpbin(:,i)]=consolidator(mdate,datp(:,i),'nanmedian',minute);  
end
datpbin(:,1)=x;

f5=figure;%overlay the binned data on the un binned
set(gcf,'OuterPosition',[100 350 1400 500]')
plot(datp(:,1),datp(:,159),'.','color',[0.7 0.7 0.7]);
hold on
line(dat(:,1),zeros(size(dat(:,1))))
title('Particulate Absorption at 678nm - minute binned')
xlabel('Julian day')
ylabel('a_{p}678 [m-1]')
hold on
plot(datp(:,1),datp(:,159),'k.');
plot(datpbin(:,1)-ref_date,datpbin(:,159),'r.') 
legend('raw-filtered (ap678)','zero (blanks should run through here)','ap678 without blanks','ap678 1-minute binned')
hold off

%% Save out
% end of this section
% saves out the important factors
cd(resultsfolder)
save(horzcat('acs_',cruiseID,'_binned','.mat'),'datpbin','ref_date','cols','cruiseID','acsfolder','band','resultsfolder');
