clear Nnew

nsbins=30;
logflag=0;

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
        x=sum(Sc(i).s)/(Send-Start)/nsbins; %number aiming for in each bin, if all are to be equal
    end
end
    
    

NnewT=0;


for i=1:length(snew)-1
 for j=1:length(Sc)
    
    Nnew(j).n(i)=0;
    suma=0;
    
        a=find( Sc(j).s>=snew(i) & Sc(j).s<snew(i+1) );
        
        la=length(a);
        if la>=1
            i1=a(1);
            i2=a(end);
            suma=sum( N(j).n(i1:i2-1) );
        else
            i1=0;
            i2=length(N(j).n)+1;
            b=find(Sc(j).s<=snew(i));
            
            if length(b)>=1 & b~=length(Sc(j).s)
                b=b(end);
                sum0=N(j).n(b) *  ( snew(i+1)-Sc(j).s(b) ) / ( Sc(j).s(b+1)-Sc(j).s(b) );
                sum02=N(j).n(b) *  ( snew(i)-Sc(j).s(b) ) / ( Sc(j).s(b+1)-Sc(j).s(b) );
                suma=sum0-sum02;
            end    
        end
        
        if i1-1>=1
            Nnew(j).n(i)= N(j).n(i1-1) - N(j).n(i1-1) *  ( snew(i)-Sc(j).s(i1-1) ) / ( Sc(j).s(i1)-Sc(j).s(i1-1) );
        end
        if i2<=length(N(j).n)
            Nnew(j).n(i) = Nnew(j).n(i) + N(j).n(i2) *  ( snew(i+1)-Sc(j).s(i2) ) / ( Sc(j).s(i2+1)-Sc(j).s(i2) );
        end
        
        Nnew(j).n(i)=Nnew(j).n(i) + suma;
        
        NnewT(i)=NnewT(i)+Nnew(j).n(i);
  end

    
end