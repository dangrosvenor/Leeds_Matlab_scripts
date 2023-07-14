topadd=26000; %height to extend sounding up to (m)
stats=[2 2]; %station numbers

pro2=pro;
for i=1:1
    
    
    
    diff=pro(2).p(2:end,1,i)-pro(2).p(1:end-1,1,i);
    izero=find(abs(diff-0)<0.01); %point where data ends in sounding
    if length(izero)==0; 
        izero=length(pro(2).p(:,1,i)) - 1;  
    else 
        izero=izero(1);
    end
    %[a closest]=min(  abs( pro(2).p(i,1,1)-pro(5).p(:,1,1) )  );
    
    
    
    izero2=find(abs(pro(stats(i)).p(:,4,i)-100)<0.01); %for 100% rh values in sounding    
    izero2=izero2(1);
	pro2(stats(i)).p(izero2:izero,10,1)=interp1(pro(5).p(:,1,i),pro(5).p(:,10,i),pro(stats(i)).p(izero2:izero,1,i),'linear','extrap');
    pro2(stats(i)).p(izero2:izero,4,1)=interp1(pro(5).p(:,1,i),pro(5).p(:,4,i),pro(stats(i)).p(izero2:izero,1,i),'linear','extrap');
    
    izero2=find(pro(stats(i)).p(:,8,i)<0.01); %points with no wind   
    pro2(stats(i)).p(izero2,8,i)=interp1(pro(5).p(:,1,i),pro(5).p(:,8,i),pro(stats(i)).p(izero2,1,i),'linear','extrap');
    pro2(stats(i)).p(izero2,9,i)=interp1(pro(5).p(:,1,i),pro(5).p(:,9,i),pro(stats(i)).p(izero2,1,i),'linear','extrap');
 
    pro2(stats(i)).p(izero+1,1,i)=topadd;
    
    for j=2:10
        pro2(stats(i)).p(izero+1,j,i)=interp1(pro(5).p(:,1,i),pro(5).p(:,j,i),topadd,'linear','extrap');
    end
    
        
       
    
end


   