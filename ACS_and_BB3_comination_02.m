close all
clear all
clc

%% ACS and BB3 combinie and finilise
% Ben Lowin
% 5/19/2021

% the goal of this script is to combine the ACS and BB3 cleaned data sets
% this is the last processing step
% maps and data ploting will be done in other scripts

%% Set up -- needs filled out
% this holds all of the varibles that controle where the data is and what
% data the program will use for this. 
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

BASE_DIR = ['F:\Back-Up(4-18-2021)\MATLAB\'];
% this should point to your MATLAB folder

Cruise_DIR=['Thesis_PhD\Cruise_Data\SKQ-June-2021\'];
% this shuold continue from the MATLAB folder to where the curise directory
% is. This should hold all of the data processed previously and the
% underway directory for this cruise.

Underway_DIR=['Underway_Data\proc\'];
% this shuold continue from the curise directory folder to where the
% underway data is. 


%% Load data sets
% This will load in the pre-made data sets
% -acs processed
% -bbp
% -topograpic map (BeringChukchiBeaufort)

cd([BASE_DIR Cruise_DIR])
%changes directory so that the acs and bb3 files are found
load (['acs_' CruiseID '_processed.mat']);
%load in the ACS file
% example- acs_SKQ_2020_jul_processed.mat
% format- acs_Cruise_ID_MMM_processed.mat

load (['bb3_' CruiseID '_processed.mat']);
% load the BB3 file
% example- bbp_SQK_2020_sep.mat
% format- bbp_CruiseID_MMM.mat

cd([BASE_DIR Cruise_DIR Underway_DIR])
%changes directory for the following data set


%% load in the underway data

UW = [];
%initilises the UnderWay matrix

%read in the underway data
a=readtable('SKQ202110S.csv','PreserveVariableNames',true);

UW=table2array(a);

% 
% for in = 1:length(days) % goes throguht each day in the underway data
%     UWF = dir([BASE_DIR Cruise_DIR Underway_DIR  days(in).name '\' cruise '_underway_' days(in).name '.txt']);
%     %opens the corosponding underway folder
%     if ~isempty(UWF) %checks if there is a matching folder
%        filename = [BASE_DIR Cruise_DIR Underway_DIR  days(in).name '\' UWF(1).name];
%        %creates a file name to evaluate the 
%        load(filename);
%        %eval(['load ' filename]);
%        %loads in the file
%        A = eval(UWF(1).name(1:end-4));  
%        %loads the data into A
%        [m,n] = size(A);
%        %determins the size of A
%       if in ==1
%           %for the first file only
%           N=n;
%           %sets the exspected file lenght
%       end
% 
%        if n == N %if the number of rows is right
%              eval(['UW = [UW;' UWF(1).name(1:end-4) '(1:end-1,:)];'])  
%              %add the new file to the matrix
%        elseif n>N
%            eval(['UW = [UW;' UWF(1).name(1:end-4) '(1:end-1,1:' num2str(N) ')];'])  
%            %if n is bigger then just take 1-N.
%            %note this causes row 37 to change... not usefull
%        else %if the number of rows is wrong - center boarrd was not deployed
%              A = [A(:,1:15) NaN*A(:,1:3) A(:,16:end)];
%              %adds in three rows for the center board
%              UW = [UW; A(1:end-1,:)];
%              %merges with underway matrix
%        end
%     end
% end
% 
% %% Quality controles the data
% % removes excess infromation
% % implments clip limits (not sure why number are what they are.
% 
% UW=UW(:,1:36);
% %removes the formast atm data
%  [m,n] = size(UW);
%  %defines the size
% 
%     for in = 1:n
%         UW(UW(:,in) < CLIPLIMITS(in,1),in) = CLIPLIMITS(in,1);
%         UW(UW(:,in) > CLIPLIMITS(in,2),in) = CLIPLIMITS(in,2);  
%         %this implments the clip limits set above 
%     end

%% Re-formats underway data into a structure
% also adds Mdate and local Mdate
%

for in=1:length(a.Properties.VariableNames)
%     NS=find(FNAME(in,:)~=' ');
%     %using the names from before find where they are not a space
    eval ([CruiseID '.' char(a.Properties.VariableNames(in)) '=UW(:,' num2str(in) ');' ]);
    %genorate the line of code to add the next layer to the structure
end


eval ([CruiseID '.mdate=datenum(UW(:,1:6));' ]);
eval ([CruiseID '.mdateloc=datenum(UW(:,1:6))-(8/24);' ]);

%% add in NaN to the time skips -ACS
% find the location of time jumps
% insert a time in there at 1 second after the end of the jump, with and
% associated NaN in the data, do this for ap, cp, lha676 and chl

diff_mdate=diff(acs.mdate);%finds the differnce between each row

diff_index=find(diff_mdate>datenum(0,0,0,0,2,0));%finds the index where more than 2 min of being off

mdate_add=acs.mdate(diff_index)+datenum(0,0,0,0,1,0);%make a time vecotr for the new data
jd_add=acs.jd(diff_index)+(1/24/60);%makes a new time vector for JD

acs_up.mdate=zeros(size(acs.mdate,1)+size(diff_index,1),size(acs.mdate,2));%creates a prealocated array
acs_up.ap=zeros(size(acs.ap,1)+size(diff_index,1),size(acs.ap,2));%creates a prealocated array
acs_up.cp=zeros(size(acs.cp,1)+size(diff_index,1),size(acs.cp,2));%creates a prealocated array
acs_up.lha676=zeros(size(acs.lha676,1)+size(diff_index,1),size(acs.lha676,2));%creates a prealocated array
acs_up.chl=zeros(size(acs.chl,1)+size(diff_index,1),size(acs.chl,2));%creates a prealocated array
acs_up.jd=zeros(size(acs.jd,1)+size(diff_index,1),size(acs.jd,2));%creates a prealocated array

acs_up.ap=nan(size(acs_up.ap));%sets the whole matrix to NaN so that the added rows are now NaN
acs_up.cp=nan(size(acs_up.cp));
acs_up.lha676=nan(size(acs_up.lha676));
acs_up.chl=nan(size(acs_up.chl));

addRows= ismember(1:size(acs.mdate,1),diff_index); %finds where the olde data goes

Olddata_index=(1:size(acs.mdate,1))+cumsum([0,addRows(1:end-1)]); %takes addes in the new data postions to that index

acs_up.mdate(Olddata_index,:)=acs.mdate(:,:);%puts the old data in
acs_up.ap(Olddata_index,:)=acs.ap(:,:);
acs_up.cp(Olddata_index,:)=acs.cp(:,:);
acs_up.lha676(Olddata_index,:)=acs.lha676(:,:);
acs_up.chl(Olddata_index,:)=acs.chl(:,:);
acs_up.jd(Olddata_index,:)=acs.jd(:,:);

newDataInd=diff_index +(1:length(diff_index))';%makes an index of the new data

acs_up.mdate(newDataInd,:)=mdate_add(:,:);%puts the new data in
acs_up.jd(newDataInd,:)=jd_add(:,:);%puts the new data in

%% match ACS to the structure 
%
% WARNING - if this code is run multiple times it will start to remove data

acs_names= strvcat('mdate','ap', 'cp',  'jd', 'lha676', 'chl', 'wavelength'); 
%the acs structure names, note the order of the names matters
%mdate must be firsts and wavelenght must be last 

for in=1:7%the ones that I want to pull out
    NS=find(acs_names(in,:)~= ' ');
    %finds the non space characters
    if in==7%does not change wavlenth
        eval([CruiseID '.acs_' acs_names(in,NS) '=acs.' acs_names(in,NS) ';'])
        %adds wavelenght unchanged to the data structure
        
    else
        eval(['acs_up.' acs_names(in,NS) '=acs_up.' acs_names(in,NS) '(1:end-1,:);' ])
        %removing the nan at the end of each data set (basicly the last
        %row)
        eval([CruiseID '.acs_' acs_names(in,NS) '=interp1(acs_up.mdate,acs_up.' acs_names(in,NS) ',' CruiseID '.mdateloc);'])
        %interplating and saving to the cruiseID structure
        
    end
end


%% add in NaN to the time skips-BB3
% find the location of time jumps
% insert a time in there at 1 second after the end of the jump, with and
% associated NaN in the data, do this for ap, cp, lha676 and chl
clear diff_mdate diff_index mdate_add jd_add addRows Olddata_index newDataInd

diff_mdate=diff(bbp.mdate);%finds the differnce between each row

diff_index=find(diff_mdate>datenum(0,0,0,0,2,0));%finds the index where more than 2 min of being off

mdate_add=bbp.mdate(diff_index)+datenum(0,0,0,0,1,0);%make a time vecotr for the new data

bbp_up.mdate=zeros(size(bbp.mdate,1)+size(diff_index,1),size(acs.mdate,2));%creates a prealocated array
bbp_up.bbp470=zeros(size(bbp.bbp470,1)+size(diff_index,1),size(bbp.bbp470,2));%creates a prealocated array
bbp_up.bbp532=zeros(size(bbp.bbp532,1)+size(diff_index,1),size(bbp.bbp532,2));%creates a prealocated array
bbp_up.bbp650=zeros(size(bbp.bbp650,1)+size(diff_index,1),size(bbp.bbp650,2));%creates a prealocated array
bbp_up.err470=zeros(size(bbp.err470,1)+size(diff_index,1),size(bbp.err470,2));%creates a prealocated array
bbp_up.err532=zeros(size(bbp.err532,1)+size(diff_index,1),size(bbp.err532,2));%creates a prealocated array
bbp_up.err650=zeros(size(bbp.err650,1)+size(diff_index,1),size(bbp.err650,2));%creates a prealocated array

bbp_up.bbp470=nan(size(bbp_up.bbp470));%sets the whole matrix to NaN so that the added rows are now NaN
bbp_up.bbp532=nan(size(bbp_up.bbp532));
bbp_up.bbp650=nan(size(bbp_up.bbp650));
bbp_up.err470=nan(size(bbp_up.err470));
bbp_up.err532=nan(size(bbp_up.err532));
bbp_up.err650=nan(size(bbp_up.err650));

addRows= ismember(1:size(bbp.mdate,1),diff_index); %finds where the olde data goes

Olddata_index=(1:size(bbp.mdate,1))+cumsum([0,addRows(1:end-1)]); %takes addes in the new data postions to that index

bbp_up.mdate(Olddata_index,:)=bbp.mdate(:,:);%puts the old data in
bbp_up.bbp470(Olddata_index,:)=bbp.bbp470(:,:);
bbp_up.bbp532(Olddata_index,:)=bbp.bbp532(:,:);
bbp_up.bbp650(Olddata_index,:)=bbp.bbp650(:,:);
bbp_up.err470(Olddata_index,:)=bbp.err470(:,:);
bbp_up.err532(Olddata_index,:)=bbp.err532(:,:);
bbp_up.err650(Olddata_index,:)=bbp.err650(:,:);

newDataInd=diff_index +(1:length(diff_index))';%makes an index of the new data

bbp_up.mdate(newDataInd,:)=mdate_add(:,:);%puts the new data in

%% BB3 data
%
%
BB3_names = strvcat('bb3_mdate', 'bb3_470', 'bb3_532', 'bb3_650', 'bb3_err470', 'bb3_err532', 'bb3_err650');
%the names to save as

BB3_struct = strvcat('mdate', 'bbp470', 'bbp532', 'bbp650', 'err470', 'err532', 'err650');
%the bbp structure names

for in=1:7
    
    NS=find(BB3_names(in,:)~=' ');
    %finds the non-space characters for the save
    NS2=find(BB3_struct(in,:)~=' ');
    %finds the non-space character for the structure
    eval([CruiseID '.' BB3_names(in, NS) '=interp1(bbp_up.mdate,bbp_up.' BB3_struct(in,NS2)  ',' CruiseID '.mdateloc);']);
    %interplating and saving to the cruiseID structure
    
end

%% Calculating Phytoplankton Carbon 
%
% 

NS=find(BB3_names(2,:)~=' ');
%makes ues of the previouse sections names to genrorate it 
eval([CruiseID '.C_phyto=12128*' CruiseID '.' BB3_names(2,NS) '+0.59;' ])
%C_phyto=12128*bbp470+0.59 --- equation taken from Graff et al 2016

%% saving the data set
%
%
%  cd ([BASE_DIR Cruise_DIR 'Results\'])
%  save([CruiseID '_full_04.mat'], [CruiseID])
