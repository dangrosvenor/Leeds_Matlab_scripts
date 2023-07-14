function data=dgrate(diag,dgs,times,dx,jjp,dgstrDan,GridDan)


dt(1)=times(1);
dt(2:length(times))=times(2:end)-times(1:end-1);



for j=1:length(diag)
    dgfind=findhead(dgs,dgstrDan(j).dg);
            
    kkp=size(diag(j).dg,1);
    dz=GridDan(j).Z(2:end)-GridDan(j).Z(1:end-1);
    dz2=repmat(dz,[1 length(dgfind)]);
    
    for t=1:length(dt)
        data(1:kkp-1,1:length(dgfind),t)=diag(j).dg(2:end,dgfind,t)*dt(t).*dz2; %multiply rate by time period and dz to get in kg/m^2
    end    
end