function [datp,yi,x159,y159]=acs_rem_blank_WB(dat,mdate,blankmdate,ref_date,blankind,cols)
j=1;
datnum=size(dat);
explim=zeros(1,2);
exp=zeros(1,2);
datp=dat;
yi=zeros(cols-1,datnum(1)); %said 165
 for p=2:cols-1
    for i=2:length(blankind)
        if blankind(i)-blankind(i-1) >100 %find start of each blank
            if p==2
                blankstart=blankmdate(i-1)-ref_date
            end
            x=blankind(j:i-1);
            y=dat(blankind(j:i-1),p);   
            imin=find(y==min(y));
            if length(imin>1)
                imin=imin(end);
            end
            if length(y)-imin>30 && imin>30 %if minimum is not near the beginning or end of the blank
                explim(1,2)=nanmedian(y(imin-30:imin+30)); %median value of the minute of data surrounding the minimum
                explim(1,1)=dat(x(imin),1); %at the minimum
            else %if minimum is near the end of the blank
                explim(1,2)=nanmedian(y(end-60:end)); %median value of last minute
                explim(1,1)=dat(x(end)-30,1); %30 seconds before end of blank 
            end
            j=i;
            exp=cat(1,exp,explim);
        end
    end
    %interpolate
    xi=mdate-ref_date;
    x=exp(2:end,1);
    y=exp(2:end,2);
    yi(p,:)=interp1(x,y,xi,'pchip','extrap');  %fit y, pchip or nearest
%     if p>85
%         hold on;
%         plot(x,y,'ro',xi,yi,'g.');
%     end
    datp(:,p)=datp(:,p)-yi(p,:)';
    j=1;
    explim=zeros(1,2);
    exp=zeros(1,2);
    if p==159
        x159=x; %so you can plot this wavelength to check if function worked
        y159=y; %so you can plot this wavelength to check if function worked
    end
end
return