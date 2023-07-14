%calcs ice MR field for all times
f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone
pp=repmat(Grid.PREFN,[1 size(potemp(1).p,2)]);

for j=1:length(potemp)
    for i=1:size(potemp(j).p,3)

        temp(j).t(:,:,i)=potemp(j).p(:,:,i)./(1e5./pp).^0.286;
        
        %T=pr(1).p(:,3)+273.15; %K
        ei=SatVapPress(temp(j).t(:,:,i),'buck2','ice'); %Pa
        %P=pr(1).p(:,2)*100; %Pa
    
        satmr(j).s(:,:,i)=f*0.622*ei./(pp-ei);
        
    end
end