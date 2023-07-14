clear Nnew NnewT Mav Mav2 macc macc2


addstartbins=1;
    
nsbins=33; %33 for 30 bins
logflag=0;
itmax=5000;

weight=0;

Sstart=min([Sc(1).s Sc(2).s]);
Send=max([Sc(1).s(end) Sc(2).s(end)]);


Send=5.5/100;



%constants to work out mass of aerosols
mu=3; %no. ions of salt that dissociate
Ms=132.1e-3; %molecular weight of salt g/mol
Mw=18.02; %molecular weight of water
T=298;
A=3.3e-7./T; %units = m
rhoS=1.769e-3; %density of salt kg/m^3 (ammonium sulphate) 
c=4/3*pi*(1/9*Ms*A^3/4.3/mu/pi);



if logflag==1 %equal in log space
Sspacing=( log(Send) - log(Sstart) ) / nsbins;
snew=[log(Sstart):Sspacing:log(Send)];
snew=exp(snew);
end

if logflag==0
	Sspacing=( (Send) - (Sstart) ) / nsbins;
	snew=[(Sstart):Sspacing:(Send)];
end

if weight==1
    x=0;
    for i=1:length(Sc)
        x=x+sum(N(i).n)/nsbins; %number aiming for in each bin, if all are to be equal
    end
end

if addstartbins==1
    weight=0;
    nex=3;
    fbin=snew30(2)-snew30(1);
    fval=snew30(1);
    snew=fval;
    for i=1:nex-1
        snew=[snew fval+i*fbin/nex];
    end
    snew=[snew snew30(2:end)];
    
end


    
    
    

%NnewT=zeros(size(N(1).n));

count=0;
it=0; %number of iterations in bisection
%for i=1:length(snew)-1
imax=length(snew)-1;
flag2=1;

nns2=[];
sss2=[];
domass=1;

while count<=imax-1
    count=count+1;
    NnewT(count)=0;
    flag=2;
    if flag2==1 %if it is the first pass for this bin
        slow=snew(count);
        shigh=snew(count+1);
    end
    flag2=1; %reset flag2
    
sss=[];
nns=[];   
  
 for j=1:length(Sc)
    
    Nnew(j).n(count)=0;
    
    suma=0;
    
        a=find( Sc(j).s>=snew(count) & Sc(j).s<snew(count+1) );
        
        la=length(a);
        if la>=1
            i1=a(1);
            i2=a(end);
            suma=sum( N(j).n(i1:i2-1) );
            
            nns=[nns N(j).n(i1:i2-1)]; %store n values onto end of vector
            sss=[sss Sc(j).s(i1:i2-1)]; %store s values
            
        else
            flag=flag-1;
            i1=0;
            i2=length(N(j).n)+1;
            b=find(Sc(j).s<=snew(count));
            
            if length(b)>=1 & b~=length(Sc(j).s)
                b=b(end);
                sum0=N(j).n(b) *  ( snew(count+1)-Sc(j).s(b) ) / ( Sc(j).s(b+1)-Sc(j).s(b) );
                sum02=N(j).n(b) *  ( snew(count)-Sc(j).s(b) ) / ( Sc(j).s(b+1)-Sc(j).s(b) );
                suma=sum0-sum02;
                
                nns=[nns suma]; %store n values onto end of vector
                sss=[sss Sc(j).s(b)]; %store s values
                
            else
                suma=0;
            end
            
            
            
        end
        
        if i1-1>=1
            Nnew(j).n(count)= N(j).n(i1-1) - N(j).n(i1-1) *  ( snew(count)-Sc(j).s(i1-1) ) / ( Sc(j).s(i1)-Sc(j).s(i1-1) );
            nns=[nns Nnew(j).n(count)];
            sss=[sss Sc(j).s(i1-1)];
        end
        if i2<=length(N(j).n)
            Nnew(j).n(count) = Nnew(j).n(count) + N(j).n(i2) *  ( snew(count+1)-Sc(j).s(i2) ) / ( Sc(j).s(i2+1)-Sc(j).s(i2) );
            nns=[nns N(j).n(i2) *  ( snew(count+1)-Sc(j).s(i2) ) / ( Sc(j).s(i2+1)-Sc(j).s(i2) )];
            sss=[sss Sc(j).s(i2)];
        end
        
        Nnew(j).n(count)=Nnew(j).n(count) + suma;
        
        NnewT(count)=NnewT(count)+Nnew(j).n(count);
  end %j loop
  
  if weight==1
        domass=0; %only calulate mass values when the bin size has been decided     
      
          
          
          
      if count+1<=length(snew)
       if NnewT(count)>x*1.1 & flag>0 & it<itmax %biection between lower and upper limits
          %imax=imax+1;
          %snew(count+1:end+1)=snew(count:end);
          shigh=snew(count+1); %refine upper limit if too high
          snew(count+1)=( shigh + slow )/2; %halve the bin
          count=count-1; %set back counter by one to retest new bin
          flag2=0;
          it=it+1;
          
         
       elseif NnewT(count)<x*0.9 & flag>0 & it<itmax
          if count==22
              '';
          end
          
           if count+2>=length(snew) 
               if abs( snew(end)-shigh ) < snew(end)/10 %if are hitting last bin
                   snew(end-1)=snew(end); 
                   count=count-1; %go through to recalculate N
                   weight=0; %but don't test size, just finish
                   snew(end)=[];
                   imax=imax-1;
                   domass=1;
               end
           elseif abs( snew(count+1)-snew(count+2) ) < snew(count+1)/100
                  shigh=snew(count+3); %if is hittinguppper boundary of bin move boundary to next bin up
                  imax=imax-1; %reduce imax
                  snew(count+1:end-1)=snew(count+2:end);%remove a bin from snew
                  snew(end)=[];
                  it=0;
                  slow=snew(count);
                  %count
           else
                 shigh=snew(count+2);
                 slow=snew(count+1);
           end
              
            
            if weight==1
              %slow=snew(count+1);
              
              snew(count+1)=( shigh + slow )/2; %halve the bin
              count=count-1; %set back counter by one to retest new bin
              flag2=0;
              it=it+1; %target cannont be reached in higher s cases due to lack of aerosol in bins
            end      
        else
            it=0; %reset it if no modifaction necessary
            
            domass=1; %work out mass later if bin width is ok
            
            
        end 
      
   end %if count+1<=length(snew) 
      
 end %weight==1
 
 if domass==1
       nns2=[nns2 nns]; %append n and s values used
       sss2=[sss2 sss];
       
       mm2=sum( c * nns ./ sss.^2 );
       mmav2(count)=mm2/sum(nns);
 end
 
  
end  %while
%snew(end)=shigh;

[sss3 isss3]=sort(sss2); %sort out all the s values used into ascending order 
nns3=nns2(isss3);  %and order n values with the s values     
massnew= c * nns3 ./ sss3.^2; %masses of all the aerosol used in above 
mnew= cumsum(massnew); %cumulative sum

mu=3; %no. ions of salt that dissociate
Ms=132.1e-3; %molecular weight of salt g/mol

Mw=18.02; %molecular weight of water

T=298;
A=3.3e-7./T; %units = m
rhoS=1.769e-3; %density of salt kg/m^3 (ammonium sulphate) 
%find dry diameters corresponding to S bin values
for i=1:length(snew)
    for j=1:length(Sc)
        if Sc(j).s(end)>=snew(i)
            [minval ind]=min(  abs( Sc(j).s-snew(i) )  ); %find s value in Sc(1).s that is closest to snew(i)
            Dnew(i)=exp(logD(j).d(ind));          
        end 
    end 
end 
 

c=4/3*pi*(1/9*Ms*A^3/4.3/mu/pi);

cc=1;
mc(1)=0;

for i=1:length(snew)-1
    %work out average mass of each bin
    Nt(i)=0;
        %estimate assuming constant dN/dD or dN/dS across bin
        if Dnew(i+1)==Dnew(i) %if d values are the same work out from constant D
            Mav(i)=pi/6*rhoS*Dnew(i)^3; 
        else
            Mav(i)=pi*rhoS/24 * (Dnew(i+1)^4 - Dnew(i)^4) / (Dnew(i+1)-Dnew(i)); %otherwise integrate assuming constant dN/dD
        end
        
        Mav2(i)=c * (1/snew(i) - 1/snew(i+1)) /(snew(i+1)-snew(i)) ; %assuming constant dN/dS
      

      for j=1:length(Sc)
          D=exp(logD(j).d);
        a=find(Sc(j).s>=snew(i) & Sc(j).s<snew(i+1));
        m=0;
        m2=0;

        for k=1:length(a)
            m(k) = c * N(j).n(a(k)) / Sc(j).s(a(k))^2;
            %m2(k) = pi*rhoS/24 * (D(a(k+1))^4 - D(a(k))^4) / (D(a(k+1))-D(a(k)));
            m2(k) = pi/6*D(a(k))^3*rhoS*N(j).n(a(k));
            cc=cc+1;
            mc(cc)=mc(cc-1)+m(k);
            sc(cc)=Sc(j).s(a(k));
        end
        mm(j)=sum(m);
        mm2(j)=sum(m2);
        Nt(i)=Nt(i)+sum(N(j).n(a));
        
      end
        
      if Nt>0
        macc(i)=sum(mm)/Nt(i);
        macc2(i)=sum(mm2)/Nt(i);
      end
        
end