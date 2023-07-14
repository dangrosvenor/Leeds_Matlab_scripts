%calcs ice MR field for all times
%f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone
%pp=repmat(Grid.PREFN,[1 502]);

casesw=6;

switch casesw
    
case 1
for j=1:length(vap)
    for i=50:size(vap(j).v,3)
        
        pcs=[0 5 10 25 50 75 90 95 100];
        
        for k=1:length(pcs)
            pcents_icemr(j).p(:,i,k)=prctile(satmr(j).s(:,:,i)',pcs(k));
        end
        
        
        
    end
end

case 2
for j=1:length(potemp)
    for i=50:size(potemp(j).p,3)
        
        pcs=[0 5 10 25 50 75 90 95 100];
        
        for k=1:length(pcs)
            pcents_potemp(j).p(:,i,k)=prctile(potemp(j).p(:,:,i)',pcs(k));
        end
        
        
        
    end
end  


case 3
    %calculate all fall speeds for ice, snow and graupel
for j=1:length(icemr)
    for i=50:size(icemr(j).i,3)
        
        TwoD.Q(:,:,6)=icemr(j).i(:,:,i);
        TwoD.Q(:,:,7)=icenc(j).i(:,:,i);
        
        Vice(j).v(:,:,i)=FallSpeed(Grid,TwoD,'ice');
        
        
    end
end      

case 4
    %calculate all fall speeds for ice, snow and graupel

for i=50:size(snowmr(1).i,3)
    for j=1:length(snowmr)
        
        TwoD.Q(:,:,4)=snowmr(j).i(:,:,i);
        TwoD.Q(:,:,9)=snownc(j).i(:,:,i);
        
        Vsnow(j).v(:,:,i)=FallSpeed(Grid,TwoD,'snow');
        
        
    end
end      
    

case 5
    
    
for j=1:length(Vsnow)
    
    b=isnan(Vsnow(j).v);
    Vsnow(j).v(b)=0;
    b=isnan(Vice(j).v);
    Vice(j).v(b)=0;


    for i=50:size(Vsnow(j).v,3)
        
        pcs=[0 5 10 25 50 75 90 95 100];
        si=size(Vsnow(j).v);
        %for k=1:length(pcs)
            pcents_vsnow(j).p(1:si(1),i,1:length(pcs))=prctile(Vsnow(j).v(:,:,i)',pcs)';
            pcents_vice(j).p(1:si(1),i,1:length(pcs))=prctile(Vice(j).v(:,:,i)',pcs)';
            %end
        
        
        
    end
end  

case 6
for j=1:length(vap)
    for i=50:size(vap(j).v,3)
        
        pcs=[0 5 10 25 50 75 90 95 100];
        
        for k=1:length(pcs)
            pcents_vap(j).p(:,i,k)=prctile(vap(j).v(:,:,i)',pcs(k));
        end
        
        
        
    end
end

    
end