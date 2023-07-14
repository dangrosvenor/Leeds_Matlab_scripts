clear diff tot

%muav=[0 cumsum(mmav2.*NnewT)];
ncum=cumsum(nns3);
for i=1:length(snew)-1
    a=find(sss3>=snew(i) & sss3<snew(i+1) );
    ncum=[0 cumsum(nns3(a(2:end)))]; %cumsum of numbers in bin i
    avcum= ( sss3(a)-sss3((a(1))) )/ ( sss3(a(end))-sss3(a(1)) ).*mmav2(i)*ncum(end); %linear mass estimate - Mav * N
    diff(i).d= abs(avcum + mnew(a(1)) - mnew(a) );
    rdiff(i).d=abs(diff(i).d/(mnew(a(1))-mnew(a))); %relative error
    
    rtot(i)=sum(rdiff(i).d)/length(rdiff(i).d);
    tot(i)=sum(diff(i).d)/length(diff(i).d); %average difference
    
    %integ(i)=trap(sss3(a),diff(i).d);
        
end

%absolute error probably more important here because are interested in total
%aerosol mass at the end - later bins contribute very little aerosol and so
%their relative error not really important