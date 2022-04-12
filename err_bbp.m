function [total_error] = err_bbp(binned_counts,err_counts,dark_counts,err_dark,scale_factor,err_scale,scat,errbeta_sw,betaP470,wall_counts,err_wall,scat_wall);
errfactor=1.1*0.04; %Uncertainty in the 1.1 is 4% (from Dall'Olmo)


%step1: C-D
err1=ones(length(binned_counts),1).*sqrt(err_counts.^2+err_dark.^2);
cminusd=binned_counts-dark_counts;
%step2: S * (C-D)
err2=sqrt((scale_factor.*err1).^2+(err_scale.*cminusd).^2);
scminusd=scale_factor.*(cminusd);
%step3: beta_sw
err3=sqrt(err2.^2+errbeta_sw(:,1).^2);
scminusdbeta=scminusd-scat; %CHECK, should be equal to betaP470

%step4: 2*pi*1.1
errtot=2*pi*sqrt((betaP470.*errfactor).^2+(1.1.*err3).^2);

%step5 - error on bbwall
cminusdwall=wall_counts-dark_counts;
err1wall=sqrt(err_wall^2+err_dark^2);
err2wall=ones(length(binned_counts),1)*sqrt((scale_factor.*err1wall).^2+(err_scale.*cminusdwall).^2);
scminusdwall=ones(length(binned_counts),1)*scale_factor.*(cminusdwall);
err3wall=sqrt(err2wall.^2+errbeta_sw(:,1).^2);
scminusdbetawall=scminusdwall-scat_wall; %equal to betaP470
errwalltot=2*pi*sqrt((scminusdbetawall.*errfactor).^2+(1.1.*err3wall).^2);

total_error=sqrt(errtot.^2+errwalltot.^2); %mean relative error = 7.8%
