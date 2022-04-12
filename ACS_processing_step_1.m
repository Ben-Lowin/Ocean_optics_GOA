close all
clear all
clc

%% ACS processing - Step 1
% 14/4/2021
% Ben Lowin

% produces a file called acs_cruiseID_compile.mat
%


%% User inputs
% You need to fill out this section with the relavant information. 
%

%1. What is your base directory for Matlab?
    base_dir='C:\Users\Benjamin Lowin\Documents\MATLAB';

%2. Where is the ACS and BB3 data you are inporting?
%the current set up, takes the base directroy and then the following info.
    acsfolder=[base_dir '\Thesis_PhD\Cruise_Data\SKQ-June-2021\ACS'];
        %this holds the .dat files
    bb3folder=[base_dir '\Thesis_PhD\Cruise_Data\SKQ-June-2021\BB3'];
        %this holds the .lvm files
%2.5 - Make a folder where you want results displayed.
    resultsfolder=[base_dir '\Thesis_PhD\Cruise_Data\SKQ-June-2021\Results'];
        
%3. What was the name, year and month of the cruise? 
    cruise='SKQ';
    year='2021'; 
    month='jun'; %month should be first three letters of the month (eg 'jan')
    
%4. What is the ACS Seriel Number? This can be found in any .dat file, but is also typically in the FILE NAME ITSELF
    % old_WB ACS=271 -- new_WB ACS=294 -- 2020_ACS=338 --
    serialnum=338;

%5. Where was the cruise?
    %Gulf of Alaska=1 --- so that code can be exspanded, currently only for
    %GOA
    location_of_cruise=1;

%% Directory and Cruise ID
% Adds all of the filles under your base directory to the path.
% Genorates and check crusie ID

%adds all folders to path
    addpath(genpath(base_dir));
    
%makes cruise ID
    cruiseID=horzcat(cruise,'_',year,'_',month);

%checks that cruise ID is correct with user.
    reply=input(['Is the following cruise ID correct? Y/N: ' cruiseID '  '] ,'s');
    
    if reply=='Y'
        disp('Excelent, moving on.')
    elseif reply=='N'
        error('Check that input information was correct.')
    else
        error('Did you use a capital letter?')
    end
    
%% Date and time
% Creates a reference data, based on the year prior to the data. 
%   eg: Dec31 of the PREVIOUS YEAR
% Finds how many hours from Grenich Mean Time the data is, based on month
%   and location.

%ref_data is the reference date in julian days
    ref_date=datenum(str2num(year)-1,12,31);
    
%loop to determin how many hours behind GMT.
%hrsbehind is the difference between where the cruise was and GMT
    if location_of_cruise==1 
        %1 mean the cruise was in the Gulf of Alaska
        if month=='jan' | month=='feb' | month=='mar' | month=='nov' ...
                | month=='dec'
            hrsbehind=9;
        elseif month=='apr'| month=='may' | month=='jun' | month=='jul' ...
                | month=='aug'| month=='sep' | month=='oct'
            hrsbehind=8;
        else
            error('Cruise month not recognized. use lower case')
        end
    else
        error(['The time zone information for cruises in that region has not been ' ...
            'set up yet.'])
        error(['Please enter the information for the time zone as compared to Grench Mean Time'])
    end

%% Check that the ACS and BB3 folder both have files in them
%This is done to prevent crashes and make sure that the data is acsssable.
%

%change directory to ACS folder
    cd(acsfolder)
%count how many .dat files there are
    nom=dir ('*.dat');
%check if data files are there and  folder is acessable
    if length(nom)<1
       error('no data files! Check your cd directory or check your \ or / in list=')
    end
%check that there are no empty data files
    for i=1:length(nom)
        if nom(i).bytes==0
            error('empty data file! Find and remove from list') %empty files will crash the acs_bin function
        end
    end

%change directory to BB3 folder
        cd(bb3folder)
%count how many .dat files there are
    nom2=dir ('*.lvm');
%check if data files are there and  folder is acessable
    if length(nom2)<1
       error('no data files! Check your cd directory or check your \ or / in list=')
    end
%check that there are no empty data files
    for i=1:length(nom2)
        if nom2(i).bytes==0
            error('empty data file! Find and remove from list') %empty files will crash the acs_bin function
        end
    end
 
    
%% sets the ACS serial number and gets the needed documentation
%As each insturment has slightly different bands.
%New codes can be added for new machines or post calibration changes.
%This also require you to make a new spreadsheet.

%code for serialnum that are not in use are in the colapsed comments. 
%{
    if str2num(year)==2016 && serialnum==271 
        hl=96;
        band= xlsread('bandname_acs271_precalib.xlsx'); %output wavelengths PRE-recalibration
    elseif serialnum==271
        hl=96;
        band= xlsread('bandname_acs271_postcalib.xlsx'); %slightly different bands for each instrument
    elseif serialnum==294
        hl=100;
        band= xlsread('bandname_acs294.xlsx'); %slightly different bands for each instrument
    end
%}    

    if serialnum==338
        hl=100;
        band= xlsread('bandname_acs338.xlsx');
    else
        error('serial number unknown')
    end

%% count number of columns in ACS data
%this is used in the ACS_bin_function
%
    columntest=textread(nom(1).name,'',-1,'headerlines',hl);
    cols=numel(columntest(1,:));    

    
    
%% Compile ACS data and bin into 1 second bins
%This function takes the mean value of every 4 rows of data.
%This effectly bins the data.
%The progress is shown in the command window (i=...)and increases with 
%every .dat file that is compiled.
%---If you see error message "out of memory", it may mean you have a .dat file
%that is empty, check filesizes are remove any files that are empty.
disp(['The program will now count up to i= ' num2str(length(nom)) '.'])

%change the directory to help with finding
    cd(acsfolder)

%runs a acs_bin_WB
    [dat] = acs_bin_WB(nom,hl,cols,ref_date);

%% Load the BB3 data
% This section will load the BB3 data, it will then normalize the number of
% columns in the data set. this is needed as the data that we collet varies
% depending on if lat, long, and other data is fed into the computer.
%The current scrip has 5 columns of data. Time, three data point and vaulve
%postion. 

%change directory to the BB3 folder
cd(bb3folder)
%list = is the information for the BB3 files
list=dir ('*.lvm');
%Initialize an empty row to build on
datbb3= zeros(1,8); 

%run through each file (ignoring the first one which seems to be a
%duplicate that is not named properly?
%does this error occour when there 
for i=2:length(list)
    %removes the header lines and saves as a temp file
    temp=textread(list(i).name,'',-1,'headerlines',22); 
    
    if length(temp(1,:))<8
        %if its shorter than 8 lines, add in the missing data as zeros
        temp=cat(2,temp,zeros(length(temp(:,1)),3));
    end
    
    if length(temp(1,:))>8 
        %if too many columns, delete all column after 8
        temp(:,9:end)=[]; 
    end

    datbb3=cat(1,datbb3,temp(3:end,:)); 
    %remove the first 2 seconds of data from each file (noisy!)
end

%These 2 lines sort the matrix by time (and remove the initial row of 0's)
[T,I] = sort(datbb3(:,1)); 
%second part removes the line of zeros
datbb3 = datbb3(I(2:end),:); 


%%
%
%
figure
set(gcf,'OuterPosition',[100 400 1400 500]')
ax(1)=subplot(2,1,1);
plot(dat(:,1),dat(:,159),'k.'); %plotting a678
xlabel('julian day')
ylabel('absorption at 678nm (a678) [m-1]')
title('Raw ACS absorption data at 678nm')
ax(2)=subplot(2,1,2);
plot(datbb3(:,1),datbb3(:,2),'b.')
hold on

plot(datbb3(:,1),datbb3(:,3),'g.')
plot(datbb3(:,1),datbb3(:,4),'r.')

xlabel('julian day')
ylabel('bb3 counts at 470nm')
title('Raw BB3 count data at 470nm')
linkaxes(ax,'x')

% the ACS and BB3 are miss alinged in time, due to the flow speed. the acs
% is the part where the blanks count. you need to know that the blanks to
% be aligned with the acs as it used the code. the bb3 tell when the vaulve
% changes. bb3 is more gradual then the acs. this needs to be done so that
% they align. this helps the QC later. 

%***Look to make sure timestamps are aligned. If not, adjust accordingly.
dat(:,1)=dat(:,1)-(65/(24*60*60)); %ACS is 65 seconds ahead of BB3???
%dat(:,1)=dat(:,1)-(5/24); %ACS is 5 hours ahead of BB3???
%keep an eye on the date and time on the computer. 
%

%%
mdate=dat(:,1)+ref_date; 
cd(resultsfolder)
save(horzcat('acs_',cruiseID,'_compile','.mat'),'dat','mdate','ref_date','cols','cruiseID','acsfolder','datbb3','band','resultsfolder');
%move the save variable after the plot - 
% should step two be just aligning the times?