close all
clear all
clc

%% BB3 processing script
% Ben Lowin
% 27/5/20212



%% Data entry
%
%

cruise = 'SKQ202110S';
% this comes from the underway files
% example- SKQ202010S_underway_20200702Z
% format- cruise_underway_dateZ
% note the letter on the end can change

CruiseID=('SKQ_2021_jun');
% this is used to find the acs and BB3 data
% example - SKQ_2020_sep'
% format - BBB_YYYY_MMM

BASE_DIR = ['C:\Users\Benjamin Lowin\Documents\MATLAB\'];
% this should point to your MATLAB folder

Cruise_DIR=['Thesis_PhD\Cruise_Data\SKQ-June-2021\'];
% this shuold continue from the MATLAB folder to where the curise directory
% is. This should hold all of the data processed previously and the
% underway directory for this cruise.

Underway_DIR=['Underway_Data\proc\'];
% this shuold continue from the curise directory folder to where the
% underway data is. 

Bad_data_made=1;
% has a bad data sheet been made for the current data set?
%if yes =1 
%if no=2

bad_data_sheet=('bad_data_SKQ_2021_jun_bb3.xlsx');
% this is the name of the bad data sheet
%if not made- comment out

%% Pathing and checking for documents
%
%

addpath(genpath([BASE_DIR]));
%adds all of the files Under the base directory to the Matlab Path

cd([BASE_DIR Cruise_DIR '\BB3']);
%changes directory to the BB3 folder
list=dir ('*.lvm');
%creates a list of the .lvm files (this is the BB3 file type
if length(list)<1 %check if there are files in the directory
    error('no data files! Check your cd directory')
end

%% Genorates the Cruise ID and ref_date
%
%

ref_date=datenum(str2num(CruiseID(5:8))-1,12,31);
%Build empty structure with this cruise ID (cruiseyear_month)


%% Underway data
% This needs to be named tsg for the BB3 procsing script
% This needs to have mdate (local time), sea surface temp, sea surface
% salinity, lat and long  (in one structure).
% note data should be minute binned as the BB3 will be matched to its time
% stamp

% Set up the underway data varibles
%{
% Column #1: Year                                                              
% Column #2: Month                                                             
% Column #3: Day                                                               
% Column #4: Hour                                                              
% Column #5: Minute                                                            
% Column #6: Second                                                            
% Column #7: ins_seapath_position - Latitude [Decimal Degrees]                 
% Column #8: ins_seapath_position - Longitude [Decimal Degrees]                
% Column #9: ins_seapath_position - Speed over ground [Knots]                  
% Column #10: ins_seapath_position - Course over Ground [Decimal Degrees True]  
% Column #11: ins_seapath_position - Heading [Degrees]                          
% Column #12: ins_seapath_position - Roll [Degrees]                             
% Column #13: ins_seapath_position - Pitch [Degrees]                            
% Column #14: ins_seapath_position - Heading [Degrees]                          
% Column #15: ins_seapath_position - Heave                                      
% Column #16: mb_em302_centerbeam - Latitude [Decimal Degrees]                  
% Column #17: mb_em302_centerbeam - Longitude [Decimal Degrees]                 
% Column #18: mb_em302_centerbeam - Depth [Meters]                              
% Column #19: met_ptu307 - Atmospheric pressure[hPa]                            
% Column #20: met_ptu307 - Air Temperature[C]                                   
% Column #21: met_ptu307 - Relative Humidity[%]                                 
% Column #22: wind_gill_fwdmast_true - Wind direction, 0 to 359 degrees,T[True] 
% Column #23: wind_gill_fwdmast_true - Wind speed, N[knots]                     
% Column #24: rad_psp-pir - LW (computed longwave downwelling irradiance,Wm**2) 
% Column #25: rad_psp-pir - SW (computed shortwave downwelling irradiance,Wm**2)
% Column #26: rad_qsr2150a - PAR [uE/m2 sec]                                    
% Column #27: flow_krohne_fwd - flow speed [m/s]                                
% Column #28: flow_krohne_fwd - volume flow [l/min]                             
% Column #29: flow_krohne_fwd - coil temperature [C]                            
% Column #30: flow_krohne_fwd - conductivity [S/m]                              
% Column #31: tsg_sbe45_fwd - Temperature [C]                                   
% Column #32: tsg_sbe45_fwd - Conductivity [S/m]                                
% Column #33: tsg_sbe45_fwd - Salinity [psu]                                    
% Column #34: tsg_sbe45_fwd - Speed of Sound [m/s]                              
% Column #35: thermo_sbe38_fwd - Intake Temperature [C]                         
% Column #36: fluoro_wetstar_fwd - Chlorophyll_a [ug/l]  
%}
for i=1
VNAME = strvcat('Year',...                                                             
'Month',...                                                             
'Day',...                                                             
'Hour',...                                                             
'Minute',...                                                             
'Second',...                                                             
'Latitude [Decimal Degrees]',...                                                             
'Longitude [Decimal Degrees]',...                                                             
'Speed over ground [Knots]',...                                                             
'Course over Ground [Decimal Degrees True]',...                                                             
'Heading [Degrees]',...                                                             
'Roll [Degrees]',...                                                             
'Pitch [Degrees]',...                                                             
'Heading [Degrees]',...                                                             
'Heave',...                                                             
'Latitude [Decimal Degrees]',...                                                             
'Longitude [Decimal Degrees]',...                                                             
'Depth [Meters]',...                                                             
'Atmospheric pressure [mbar]',...                                                             
'Air Temperature [^\circC] ',...                                                             
'Relative Humidity [%] ',...                                                             
'Wind direction [^\circTrue]',...                                                             
'Wind speed [kt] ',...                                                             
'Longwave downwelling irradiance [W m^-^2]',...                                                             
'Shortwave downwelling irradiance [W m^-^2]',...                                                             
'PAR [uE m^-^2 s^-^1]',...                                                             
'flow speed [m s^-^1]',...                                                             
'Seawater Flow [l min^-^1] ',...                                                             
'coil temperature [^\circC]',...                                                             
'conductivity [S/m] ',...                                                             
'Temperature [^\circC] ',...                                                             
'Conductivity [S/m]',...                                                             
'Salinity [psu]  ',...                                                             
'Speed of Sound [m/s]',...                                                             
'Seachest Intake Temperature [^\circC]',...
'Chlorophyll_a [ug/l] ',... 
'Surface Nitrate [uM]'); 

FNAME = strvcat('Year',...                                                             
'Month',...                                                             
'Day',...                                                             
'Hour',...                                                             
'Minute',...                                                             
'Second',...                                                             
'Latitude',...                                                             
'Longitude',...                                                             
'SOG',...                                                             
'COG',...                                                             
'Heading',...                                                             
'Roll',...                                                             
'Pitch',...                                                             
'Heading2',...                                                             
'Heave',...                                                             
'Latitude',...                                                             
'Longitude',...                                                             
'Depth',...                                                             
'SLP',...                                                             
'AirTemperature',...                                                             
'RH',...                                                             
'WindDir',...                                                             
'WindSpd',...                                                             
'DLWR',...                                                             
'SDWR',...                                                             
'PAR',...                                                             
'SeachestFlow',...                                                             
'SeachestFlow2',...                                                             
'CoilTemp',...                                                             
'SeacestCond',...                                                             
'Temperature',...                                                             
'Conductivity',...                                                             
'Salinity',...                                                             
'SoundSpd',...                                                             
'SeachestIntakeTemp',...
'Chla',...
'NO3');                         

CLIPLIMITS = [1900 2100;...                                                             
1 12;...                                                             
1 30;...                                                             
0 23;...                                                             
0 59;...                                                             
0 60;...                                                             
0 90;...                                                             
-360 360;...                                                             
0 20;...                                                             
0 360;...                                                             
0 360;...                                                             
-360 360;...                                                             
-360 360;...                                                             
0 360;...                                                             
0 100;...                                                             
0 90;...                                                             
-360 360;...                                                             
0 1000000;...                                                             
900 1100;...                                                             
-100 50;...                                                             
0 150;...     % RH                                                        
0 360;...                                                             
0 150;...                                                             
-10 4000;...                                                             
-10 4000;...                                                             
0 10000;...                                                             
0 10;...                                                             
0 10;...                                                             
-2 18 ;...                                                             
0 5;...                                                             
-2 18;...                                                             
0 5;...                                                             
10 40;...                                                             
1000 2000;...                                                             
-2 20;...
-1 100;...
-10 25];                         
end %this is a loop to hid the huge number of lines used so that it is read able

cd([BASE_DIR Cruise_DIR Underway_DIR])

% load in the underway data
days = dir([BASE_DIR Cruise_DIR Underway_DIR cruise(4:7) '*']);
% days is a structure with 6 feilds. 
%{
% name
% folder
% date
% bytes
% isdir
% datenum
%}
UW = [];
%initilises the UnderWay matrix

for in = 1:length(days) % goes throguht each day in the underway data
    UWF = dir([BASE_DIR Cruise_DIR Underway_DIR  days(in).name '\' cruise '_underway_' days(in).name '.txt']);
    %opens the corosponding underway folder
    if ~isempty(UWF) %checks if there is a matching folder
       filename = [BASE_DIR Cruise_DIR Underway_DIR  days(in).name '\' UWF(1).name];
       %creates a file name to evaluate the 
       load(filename)
%      eval(['load (' filename ')']);
       %loads in the file
       A = eval(UWF(1).name(1:end-4));  
       %loads the data into A
       [m,n] = size(A);
       %determins the size of A
      if in ==1
          %for the first file only
          N=n;
          %sets the exspected file lenght
      end

       if n == N %if the number of rows is right
             eval(['UW = [UW;' UWF(1).name(1:end-4) '(1:end-1,:)];'])  
             %add the new file to the matrix
       elseif n>N
           eval(['UW = [UW;' UWF(1).name(1:end-4) '(1:end-1,1:' num2str(N) ')];'])  
           %if n is bigger then just take 1-N.
           %note this causes row 37 to change... not usefull
           
       else %if the number of rows is wrong - center boarrd was not deployed
             A = [A(:,1:15) NaN*A(:,1:3) A(:,16:end)];
             %adds in three rows for the center board
             UW = [UW; A(1:end-1,:)];
             %merges with underway matrix
       end
    end
end

% Quality controles the data
% removes excess infromation
% implments clip limits (not sure why number are what they are.

UW=UW(:,1:36);
%removes the formast atm data
 [m,n] = size(UW);
 %defines the size

    for in = 1:n
        UW(UW(:,in) < CLIPLIMITS(in,1),in) = CLIPLIMITS(in,1);
        UW(UW(:,in) > CLIPLIMITS(in,2),in) = CLIPLIMITS(in,2);  
        %this implments the clip limits set above 
        
    end

% Re-formats underway data into a structure
% also adds Mdate and local Mdate
%

for in=[1:8,31,33]
    NS=find(FNAME(in,:)~=' ');
    %using the names from before find where they are not a space
    eval (['tsg_p.' FNAME(in,NS) '=UW(:,' num2str(in) ');' ]);
    %genorate the line of code to add the next layer to the structure
end
eval (['tsg_p.mdate=datenum(UW(:,1:6));' ]);
eval (['tsg_p.mdateloc=datenum(UW(:,1:6))-(8/24);' ]);

% rename varables so script can read them

tsg.mdate=tsg_p.mdateloc;
tsg.sal=tsg_p.Salinity;
tsg.sst=tsg_p.Temperature;
tsg.lat=tsg_p.Latitude;
tsg.long=tsg_p.Longitude;


%% Time difference between computer and underway system
% This is assumed to be 0 seconds
%
% -'ve = TSG clock AHEAD of BB3 computer (e.g. TSG in UTC, BB3 in PST), so SUBTRACT time to match BB3 time
% +'ve = TSG clock BEHIND BB3, so need to ADD time to get to BB3 time
timediff=-datenum(0,0,0,0,0,0); 

x=input(['Is the following time difference correct ' num2str(timediff) ' ? If Yes enter Y, else enter N     ' ], 's');

if x=='N'
    error('time difference has not be set correctly')
elseif x=='Y'
    disp('continuing')
end

%% Filter change + lense + flags
%
%
x=input('Have you changed madatory flags? If Yes enter Y, else enter N   ' ,'s');
if x=='N'
    error('flags needs to be changed.')
elseif x=='Y'
    disp('continuing')
end


filt.mdate(1)=datenum(2021,06,27,12,00,00);

lens.mdate(1)=datenum(2021,06,27,12,00,00);

%NOTE: At least 1 flag is needed, so use beginning of cruise as 1st default 
flag.start(1)=datenum(2021,06,27,12,00,00);
flag.end(1)=datenum(2021,06,27,12,00,00);
flag.reason(1)=cellstr('Start of Cruise');

%% Plug in instrument/casket/cruise specific information:
%Which BB3 did you use? old BB3=1415 new BB3=1615
serialnum=6077; %*** CHANGE THIS

if serialnum==6077
    %Uncertainty in counts (C)
    inst.errcnts470=1; %Stated WETLABS 'Instrument resolution'
    inst.errcnts532=1.6;
    inst.errcnts650=1.5;
    %Scaling Factors (constants for your BB3)
    inst.SF470=0.00001148; 
    inst.SF532=0.000007978; 
    inst.SF650=0.000004051; 
    
elseif serialnum==1615
    %Uncertainty in counts (C)
    inst.errcnts470=1.0; %Stated WETLABS 'Instrument resolution'
    inst.errcnts532=1.2;
    inst.errcnts650=1.4;
    %Scaling Factors (constants for your BB3)
    inst.SF470=0.00001089; 
    inst.SF532=0.000007057; 
    inst.SF650=0.000003587; 
end

    inst.errSF470=1.2e-7; %from Dall'Olmo
    inst.errSF532=0.31e-7;%from Dall'Olmo
    inst.errSF650=inst.SF650*mean([inst.errSF470/inst.SF470,inst.errSF532/inst.SF532]); %SF650 *  average of relative errors for 470 and 532;

x=input('Have you changed the dark count data? If Yes enter Y, else enter N   ' ,'s');
if x=='N'
    error('dark count needs to be changed.')
elseif x=='Y'
    disp('continuing')
end
    
%Dark Counts (in casket, w/ seawater, measured during cruise):
darkcnts.cnts470=50; %xxx dark counts throughout cruise
darkcnts.cnts532=45; %xxx
darkcnts.cnts650=46; %xxx 
%Uncertainty in dark counts (D)
darkcnts.err470=1; % 0.6=stdev of 3 dark count measurements during cruise
darkcnts.err532=5; % 1.4 See DarkCounts_lip0217 for details
darkcnts.err650=2; % 1.9

%Wall Counts (measured in the lab prior to cruise) 
%**
wall.wall470=73; %
wall.wall532=77;
wall.wall650=79;
wall.err470=2;
wall.err532=2;
wall.err650=2;
wall.sal=0;
wall.temp=23;

%%
%
%
cd([BASE_DIR Cruise_DIR '\BB3\']);
[bbp]=bbp_processor_complete(list,ref_date,tsg,cd,timediff,filt,lens,flag,inst,darkcnts,wall);

%% Quality control
%
%
figure (5)
plot(bbp.mdate-ref_date,bbp.bbp470)
hold on
plot(bbp.mdate-ref_date,bbp.bbp532)
plot(bbp.mdate-ref_date,bbp.bbp650)
hold off


if Bad_data_made==1
    [num,txt,~] = xlsread(bad_data_sheet);
    for i=1:length(num(:,1))
        baddata.start(i)=num(i,1);
        baddata.end(i)=num(i,2);
        baddata.reason(i,:)=txt(i+1,3)'; 
        inmdate=find(bbp.mdate-ref_date> baddata.start(i) & bbp.mdate-ref_date<baddata.end(i));
 %remove from data
        bbp.bbp470(inmdate)=nan;
        bbp.bbp532(inmdate)=nan;
        bbp.bbp650(inmdate)=nan;
        bbp.err470(inmdate)=nan;
        bbp.err532(inmdate)=nan;
        bbp.err650(inmdate)=nan;
    end
else
    error('Pause here and make bad_data.xlsx sheet')
end

figure (6)
plot(bbp.mdate-ref_date,bbp.bbp470)
hold on
plot(bbp.mdate-ref_date,bbp.bbp532)
plot(bbp.mdate-ref_date,bbp.bbp650)
hold off



%% save out the bbp data
%
%

cd([BASE_DIR Cruise_DIR '\Results\']);


save( horzcat('bb3_',CruiseID,'_processed','.mat'),'bbp','ref_date');


