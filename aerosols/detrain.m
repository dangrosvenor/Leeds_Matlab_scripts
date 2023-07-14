function [data,flux]=detrain(diag,dgs,times,dgstrDan)


a=find(times~=0);

dt(2:length(times))=times(2:end)-times(1:end-1);
dt(a(1))=dt(a(2));

for j=1:length(diag)
    
    kkp=size(diag(j).dg,1);
    
dgfind=findhead(dgs,dgstrDan(j).dg);

for k=1:length(dgfind)
    
    
    dg_str=dgstrDan(j).dg(dgfind(k));
    
    
    if length(dg_str{1})>=4
        if strcmp(dg_str{1}(4),'_')==1
            area_str=strcat(dg_str{1}(1:3),'_A');
            dgfind2=findhead(area_str,dgstrDan(j).dg);
        end
    else
        dg_str='xxxx';
    end
    
    
    for t=1:size(diag(j).dg,3);
        
        if strcmp(dg_str{1}(4),'_')==1
            
            area(1:kkp,k)=diag(j).dg(:,dgfind2(1),t);
            areazeroes=find(area(:,k)==0);
            area(areazeroes,k)=1;

            
        else
            area(1:kkp,k)=ones([kkp 1]);
        end
    
    
        flux(1:kkp,k,t)= diag(j).dg(1:end,dgfind(k),t) * dt(t) ;%./ area(1:kkp,k);
        data(1:kkp-1,k,t)=( -diag(j).dg(2:end,dgfind(k),t)+diag(j).dg(1:end-1,dgfind(k),t) )*dt(t) ;%./ area(2:kkp,k); %mass of material detrained at this height
    
    end
end


end
