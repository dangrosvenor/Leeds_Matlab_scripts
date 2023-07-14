clear diff c
            tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]); %ref potemp
            T=squeeze(TwoD(idir).TH1)+tref; %potemp
            pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
            P=TwoD(idir).PP; %actual p
            T=T./(1e5./P).^0.286; %actual T
            tref=tref./(1e5./pref).^0.286; %ref temp
            
            Tp=T-tref; %perturbation of temperature
            

ih1=findheight(GridDan(idir).Z+620,15e3);
ih2=findheight(GridDan(idir).Z+620,17e3);

    
idir=1;
b=[];
for i=1:size(GridDan(idir).Y1)
    sig=sign(Tp(ih1:ih2,i));
    dif=diff(sig);
    [a]=find(abs(dif)==2); %where changes sign
    if length(a)>0
        b=[b GridDan(idir).Z(min(a)+ih1)+620]; %add new height onto end of b array
        c(i)=GridDan(idir).Z(min(a)+ih1)+620;
    else
        c(i)=NaN;
    end
    
end

h=mean(b)
    
    