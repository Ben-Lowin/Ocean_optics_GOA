function [bbp] = bbp_processor_complete(list,ref_date,tsg,cd,timediff,filt,lens,flag,inst,darkcnts,wall);
%% (1) Compile raw data, sort by time, and plot                   

datbb3= zeros(1,8); %Initialize an empty row to build on

for i=1:length(list)
    temp=textread(list(i).name,'',-1,'headerlines',22); 
    if length(temp(1,:))<8
        temp=cat(2,temp,zeros(length(temp(:,1)),3));
    end
    if length(temp(1,:))>8 %if too many columns, delete all column after 8
        temp(:,9:end)=[]; 
    end
    
    datbb3=cat(1,datbb3,temp(3:end,:)); %remove the first 2 seconds of data from each file (noisy!)
end

[T,I] = sort(datbb3(:,1)); %These 2 lines sort the matrix by time (and remove the initial row of 0's)
datbb3 = datbb3(I(2:end),:);

for i=1:length(datbb3(:,1))-1
    time_interval_bb3(i)= datbb3(i+1,1)-datbb3(i,1); %consider removing this
end

clearvars T I i I list temp time_interval_bb3

figure
plot(datbb3(:,1),datbb3(:,2),'b.-'); %Plot the raw data
xlabel('Time (m-date)'); 
ylabel('raw counts');
title('Figure 0: Raw counts at 470nm')
hold on
plot(datbb3(:,1),datbb3(:,3),'g.-');
plot(datbb3(:,1),datbb3(:,4),'r.-');

figure
plot(datbb3(:,1)+ref_date,datbb3(:,2),'k.-'); %Plot the raw data

xlabel('Time (dd-mm HH:MM local time)'); 
ylabel('raw counts');
title('Figure 1: Raw counts at 470nm')
datetick('x','dd-mm HH:MM','keeplimits')
z = zoom(gcf); p = pan(gcf); set(z,'ActionPostCallback',@zoomDateTick); set(p,'ActionPostCallback',@panDateTick); 

% Assign Variables
time=datbb3(:,1); %this is in Julian Day (local time) 

cnts470=datbb3(:,2);
cnts532=datbb3(:,3);
cnts650=datbb3(:,4);
valve=datbb3(:,5);
mdate=time+ref_date;
round(datevec(mdate)); %round data to every second, so they don't clump together
mdate=datenum(ans);
cnts470raw=cnts470;
cnts532raw=cnts532;
cnts650raw=cnts650;
mdateraw=mdate;
%% (2) Find blank periods, plot them, and then remove them
blanks=find(valve==1);
blanktime=mdate(blanks);
blanks470=cnts470(blanks);
blanks532=cnts532(blanks);
blanks650=cnts650(blanks);

%Overlay blanks on Figure 1
hold on
plot(blanktime,blanks470,'g.')

lensclean = 78*ones(length(lens.mdate),1); %just to plot them in a good place
plot(lens.mdate,lensclean,'m*','MarkerSize',10);

filterchange = 75*ones(length(filt.mdate),1); %just to plot them in a good place
plot(filt.mdate,filterchange,'m^','MarkerSize',10);

for i=1:length(flag.start) %Plot Flags to help locate areas of bad data
    plot([flag.start(i) flag.end(i)],75*ones(2),'r-','linewidth',2);
    text(flag.start(i),70,flag.reason{i},'fontsize',10,'color','r')

end

legend('raw data','filter periods (blanks)','lens cleaned','filter changed','flags')


% REMOVE BLANKS AND ZERO-COUNTS 

cnts470(blanks)=nan;
cnts532(blanks)=nan;
cnts650(blanks)=nan;

%Remove data 20 seconds prior to blank starting (10 seconds of transition +
%10 seconds of odd data prior to transition.
%AND 150 seconds after blank ended (chamber must refill with raw sw)
p2 = find(valve==2);   
for pp = 1:length(p2)
    cnts470(p2(pp)-10:p2(pp)+150) = nan;
    cnts532(p2(pp)-10:p2(pp)+150) = nan;
    cnts650(p2(pp)-10:p2(pp)+150) = nan;
end
%In case above statement artificially lengthens cnts vectors
cnts470=cnts470(1:length(mdate));
cnts532=cnts532(1:length(mdate));
cnts650=cnts650(1:length(mdate));

% Remove zero-counts (e.g. IF LABVIEW is running while instrument is unplugged)
cnts470(find(cnts470==0))=nan;
cnts532(find(cnts532==0))=nan;
cnts650(find(cnts650==0))=nan; 

%% (3) Bin data, incorporate TSG data, and plot


minute=1/24/60; %CHANGE HERE for different size bins
[mdatebin,cnts470bin]= consolidator(mdate,cnts470,'nanmedian',minute); 
[mdatebin,cnts532bin]= consolidator(mdate,cnts532,'nanmedian',minute);
[mdatebin,cnts650bin]= consolidator(mdate,cnts650,'nanmedian',minute);

%make artificial TSG data if TSG is not available
if length(tsg.mdate)==1
    tsg.mdate=mdatebin;
    tsg.sal=ones(length(tsg.mdate),1)*32;
    tsg.sst=ones(length(tsg.mdate),1)*15;
    tsg.lat=ones(length(tsg.mdate),1)*49;
    tsg.long=ones(length(tsg.mdate),1)*(-130);
end

tsg.mdatematch=tsg.mdate+timediff; %this makes sure tsg.mdate and mdate are the same

lat=interp1(tsg.mdatematch,tsg.lat,mdatebin);
long=interp1(tsg.mdatematch,tsg.long,mdatebin);
sst=interp1(tsg.mdatematch,tsg.sst,mdatebin);
sal=interp1(tsg.mdatematch,tsg.sal,mdatebin);


%
figure
plot(mdateraw,cnts470,'k.-');
hold on
plot(mdatebin,cnts470bin,'r.','markersize',6); %check that the +/- sign of timediff is correct here
set(gca,'FontSize',10);
set(gca,'FontSize',10);
title('Figure 2: Raw vs. Binned Counts at 470nm')
xlabel('Time [mm-dd HH:MM local time]','fontweight','bold');
ylabel('counts');
legend('Raw 470nm Counts','1-Min Median Binned')
datetick('x','dd-mm HH:MM','keeplimits','keepticks')
z = zoom(gcf);
p = pan(gcf);
set(z,'ActionPostCallback',@zoomDateTick);
set(p,'ActionPostCallback',@panDateTick);

%% (4) Convert to beta, compute wall effect, calculate bbp, and plot
%First, convert to beta using cruise-specific dark counts and intrument-specified scaling factor (SF)
beta470=(cnts470bin-darkcnts.cnts470)*inst.SF470;
beta532=(cnts532bin-darkcnts.cnts532)*inst.SF532;
beta650=(cnts650bin-darkcnts.cnts650)*inst.SF650;


wavelengths=[470 532 650];
for p=1:length(wavelengths)
    
for i=1:length(sal)
  [theta,betasw,bsw,beta90sw]=betasw_ZHH2009(wavelengths(p),sal(i),sst(i));
  scat(i,p)=betasw(11707);
  errbeta_sw(i,p)=scat(i,p).*0.0224; %error is 2.24% from Zhang et al. (2009)
end
end

betaP470=beta470-scat(:,1);
betaP532=beta532-scat(:,2);
betaP650=beta650-scat(:,3);

% Compute the scattering by the chamber walls (bbwall)
%First, use the wall counts to compute beta_wall
betawall470=(wall.wall470-darkcnts.cnts470)*inst.SF470;
betawall532=(wall.wall532-darkcnts.cnts532)*inst.SF532;
betawall650=(wall.wall650-darkcnts.cnts650)*inst.SF650;
%Compute the pure water component of the wall scattering
wavelengths=[470 532 650];
for p=1:length(wavelengths)
  
  [theta,betasw,bsw,beta90sw]=betasw_ZHH2009(wavelengths(p),wall.sal,wall.temp,0.039);
  scatwall(1,p)=betasw(11707);
end
%compute the particulate backscatter by the wall
bbwall470=(2*pi*1.1)*(betawall470-scatwall(1));
bbwall532=(2*pi*1.1)*(betawall532-scatwall(2));
bbwall650=(2*pi*1.1)*(betawall650-scatwall(3));

%compute bbp by subtracting bb_wall
bbp470=((2*pi*1.1)*betaP470)-bbwall470;
bbp532=((2*pi*1.1)*betaP532)-bbwall532;
bbp650=((2*pi*1.1)*betaP650)-bbwall650;

%
figure
plot(mdatebin,bbp470,'b.');
hold on
plot(mdatebin,bbp532,'g.');
plot(mdatebin,bbp650,'r.');
xlabel('time (dd-mm HH:MM local time)')
ylabel('bbp [m^{-1}]');
set(gca,'FontSize',10)
legend('470nm','532nm','650nm','Location','NorthWest');
datetick('x','dd-mm HH:MM','keeplimits','keepticks')
z = zoom(gcf); p = pan(gcf); 
set(z,'ActionPostCallback',@zoomDateTick); 
set(p,'ActionPostCallback',@panDateTick);

%% (5) Calculate Total Uncertainty in bbp
%Equation 1 (Dall'Olmo, 2009)
%bbp=2*pi*1.1*(S(C-D)-beta_sw)-bb_wall

%run total error function.
err470= err_bbp(cnts470bin,inst.errcnts470,darkcnts.cnts470,darkcnts.err470,inst.SF470,inst.errSF470,scat(:,1),errbeta_sw(:,1),betaP470,wall.wall470,wall.err470,scatwall(1));
err532= err_bbp(cnts532bin,inst.errcnts532,darkcnts.cnts532,darkcnts.err532,inst.SF532,inst.errSF532,scat(:,2),errbeta_sw(:,2),betaP532,wall.wall532,wall.err532,scatwall(2));
err650= err_bbp(cnts650bin,inst.errcnts650,darkcnts.cnts650,darkcnts.err650,inst.SF650,inst.errSF650,scat(:,3),errbeta_sw(:,3),betaP650,wall.wall650,wall.err650,scatwall(3));

%clearvars -EXCEPT 
    %betaP650 beta90sw beta470_tsg beta532_tsg beta650_tsg theta
    
bbp.bbp470=bbp470;
bbp.bbp532=bbp532;
bbp.bbp650=bbp650;
bbp.err470=err470;
bbp.err532=err532;
bbp.err650=err650;
bbp.cnts470=cnts470raw;
bbp.cnts532=cnts532raw;
bbp.cnts650=cnts650raw;
bbp.valve=valve;
bbp.mdateraw=mdate;
bbp.mdate=mdatebin;
bbp.lat=lat;
bbp.long=long;
bbp.sal=sal;
bbp.sst=sst;
bbp.blankmdate=blanktime;
bbp.blankcnts470=blanks470;
end
