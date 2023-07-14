%calcs ice MR field for all times
f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone
%pp=repmat(Grid.PREFN,[1 size(potemp,2)]);

for j=1:length(potemp)
    for i=44:size(potemp(j).p,3)

        temp(j).t(:,:,i)=potemp(j).p(:,:,i)./(1e5./pressure(j).p(:,:,i)).^0.286;
        ei=SatVapPress(temp(j).t(:,:,i),'buck2','ice'); %Pa
        satmr(j).s(:,:,i)=f*0.622*ei./(pressure(j).p(:,:,i)-ei);
        
    end
end

%         j=1;

%         temp(j).t=potemp(j).p./(1e5./pressure(j).p).^0.286;
%         ei=SatVapPress(temp(j).t,'buck2','ice'); %Pa
%         satmr(j).s=f*0.622*ei./(pressure(j).p-ei);