%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acs bin into 1min
% by Will Burt
% 11th, Jul, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dat] = acs_bin_WB(nom,hl,cols,ref_date)
    %read name
    dat= zeros(1,cols); %initialize dat using column information
    for i=1:length(nom)
        i
        afid = fopen(nom(i).name,'rt'); %open .dat file
        %scan text file: there are 167 columns of data in acs file (including time stamp, abs. and atten. data)
        %'headerlines' represents the number of rows chopped off at the
        %top ... for the acs data, this is always 97
        acs = textscan(afid,repmat('%f',1,cols),1000000000,'headerlines',hl,'collectoutput',1); 
        acs = acs{1};
        fclose(afid);
        timestamp=acs(:,1); %original timestamp in weird millisecond format

        name=nom(i).name;
        yyyy=str2num(name(9:12));
        month=str2num(name(13:14));
        dd=str2num(name(15:16));
        hh=str2num(name(17:18));
        mm=str2num(name(19:20));
        ss=str2num(name(21:22));
        date=datenum(yyyy,month,dd)-ref_date;
        time=(hh*3600+mm*60+ss)/86400;
        filestart=date+time;
        jd=zeros(length(timestamp)-1,1);
        for j=1:length(timestamp)
        jd(j)=filestart+(timestamp(j)-timestamp(1))/86400000;
        end
       
      second=1/24/60/60;
        [x,y]=consolidator(jd,acs(:,2),'nanmean',second);  %find out how many rows you'll need
        acsbin=zeros(length(x),length(acs(1,:)));
        
      for k=2:length(acs(1,:))-1
      [x,acsbin(:,k)]=consolidator(jd,acs(:,k),'nanmean',second);  
      end
      acsbin(:,1)=x;  
        dat=cat(1,dat,acsbin);
    end
    dat=dat(2:end,:); %remove the initial line of zeros
    [t,i]=sort(dat);
    dat=dat(i(:,1),:);
return
