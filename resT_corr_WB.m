function [ap_cor,cp_cor]=resT_corr_WB(dat_a,dat_c,Fai_S)
    %%Residual_T correction
    

    %ref is [156:171], but [170,171] is Nan
    %NIR is [166]
    b_ref= dat_c(:,156:169)-dat_a(:,156:169);
    ap_ref=dat_a(:,156:169);
    aw_ref=Fai_S(156:169,2); %according to their codes, they used 2nd colunm
    ap_nir=dat_a(:,166)*ones(1,14);
    aw_nir=Fai_S(166,6)*ones(1,14);
    b_nir=(dat_c(:,164)-dat_a(:,164))*ones(1,14);

    %for ap correction
    aw_n = Fai_S(:,2);
    b_n = dat_c-dat_a;
    ap_nir_n=dat_a(:,166)*ones(1,176);
    aw_nir_n=Fai_S(166,6)*ones(1,176);
    b_nir_n=(dat_c(:,164)-dat_a(:,164))*ones(1,176);

    %GET minimal
    dT=zeros(size(dat_a,1),1);
    ap_cor = zeros(size(dat_a,1),176);
    for i=1:length(b_ref)
        %x0=0;
        i
        options=optimset('Display', 'off', 'MaxIter' , 20000000, 'MaxFunEvals', 20000, 'TolX', 1e-8, 'TolFun', 1e-8 );
        [dT(i,1), fval, exitflag]=fminsearch(@(x)sum(abs(ap_ref(i,:)-aw_ref'.*x-(ap_nir(i,:)-aw_nir.*x)./b_nir(i,:).*b_ref(i,:))),...
        [0],options);

        %correction
        ap_cor(i,:)=dat_a(i,:)-aw_n'.*dT(i,1)-(ap_nir_n(i,:)-aw_nir_n.*dT(i,1))./b_nir_n(i,:).*b_n(i,:);
    end
    cp_cor=ap_cor+b_n;
    %ap_cor is the absorption data after all correction!
return