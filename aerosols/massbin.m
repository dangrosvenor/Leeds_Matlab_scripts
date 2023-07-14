Combined_Sdist4;

tol=6.5e-17;

tol=2e-17; %2.5e-17
weight=0;
massbinflag=1;
j=0;

while massbinflag==1
    j=j+1;
    masserr;
    iab=find(tot>tol);
    if length(iab)<1;massbinflag=0;break;end
    
    for i2b=1:length(iab)
        i=iab(i2b);
        snew(i+2:end+1)=snew(i+1:end); %make space for new bin
        snew(i+1)=(snew(i)+snew(i+1)) /2; %new bin halfway
        iab(i2b+1:end)=iab(i2b+1:end)+1;
    end
    
    Combined_Sdist_choose; %test new snew bin values
    
end

for i=1:length(snew)
    [a b]=min(abs(Sc(1).s-snew(i)));
    [a2 b2]=min(abs(Sc(2).s-snew(i)));
    if a<a2
        Ds(i)=exp(interp1(Sc(1).s,logD(1).d,snew(i)));
    else
        Ds(i)=exp(interp1(Sc(2).s,logD(2).d,snew(i)));
    end
end

Dsmid=(Ds(2:end)+Ds(1:end-1))/2;


rhoS=1.769e-3; %density of salt kg/m^3 (ammonium sulphate) 

Ms=1/6*pi*rhoS*Ds.^3;