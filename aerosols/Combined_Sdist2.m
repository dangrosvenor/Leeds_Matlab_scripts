clear Nnew NnewT

nsbins=21;
logflag=0;

weight=0;

Sstart=min([s1(1) s2(1)]);
%Send=max([s1(end) s2(end)]);

%Send=s1(end);
Send=5.5/100;

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
    
    

%NnewT=zeros(size(N(1).n));

count=0;
%for i=1:length(snew)-1
imax=length(snew)-1;

while count<=imax-1
    count=count+1;
    NnewT(count)=0;
    flag=2;
    
 for j=1:length(Sc)
    
    Nnew(j).n(count)=0;
    
    suma=0;
    
        a=find( Sc(j).s>=snew(count) & Sc(j).s<snew(count+1) );
        
        la=length(a);
        if la>=1
            i1=a(1);
            i2=a(end);
            suma=sum( N(j).n(i1:i2-1) );
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
            end    
        end
        
        if i1-1>=1
            Nnew(j).n(count)= N(j).n(i1-1) - N(j).n(i1-1) *  ( snew(count)-Sc(j).s(i1-1) ) / ( Sc(j).s(i1)-Sc(j).s(i1-1) );
        end
        if i2<=length(N(j).n)
            Nnew(j).n(count) = Nnew(j).n(count) + N(j).n(i2) *  ( snew(count+1)-Sc(j).s(i2) ) / ( Sc(j).s(i2+1)-Sc(j).s(i2) );
        end
        
        Nnew(j).n(count)=Nnew(j).n(count) + suma;
        
        NnewT(count)=NnewT(count)+Nnew(j).n(count);
  end
  
  if weight==1  
      if NnewT(count)>x & flag>0
          imax=imax+1;
          snew(count+1:end+1)=snew(count:end);
          snew(count+1)=( snew(count+2) + snew(count) )/2; %halve the bin
          count=count-1; %set back counter by one to retest new bin
          
      end
  end
  
end %while